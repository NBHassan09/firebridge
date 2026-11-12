"""Shared code for near-real-time (NRT) VIIRS hotspots from NASA GIBS vector tiles (no API key).
Used once to backfill 2026 and daily by the GitHub Action."""
import gzip, io, json, math, struct, urllib.request
import numpy as np

GIBS = "https://gibs.earthdata.nasa.gov/wmts/epsg4326/best/VIIRS_{sat}_Thermal_Anomalies_375m_All/default/{day}/500m/{z}/{r}/{c}.mvt"
Z, TILES = 6, [(14, 59), (14, 60), (15, 59), (15, 60)]          # 4.5-degree tiles covering Bangladesh
BBOX = (20.3, 26.8, 87.9, 92.8)                                    # lat_min, lat_max, lon_min, lon_max
SAT_CODE = {"SNPP": "VS", "NOAA20": "V1", "NOAA21": "V2"}

# ---------- Mapbox vector tile decoding (points only) ----------
def _varint(b, p):
    r = s = 0
    while True:
        x = b[p]; p += 1; r |= (x & 0x7F) << s; s += 7
        if not x & 0x80: return r, p
def _msg(b, p, end):
    while p < end:
        tag, p = _varint(b, p); f, w = tag >> 3, tag & 7
        if w == 0: v, p = _varint(b, p); yield f, w, v, p
        elif w == 2:
            n, p = _varint(b, p); yield f, w, (p, p + n), p + n; p += n
        elif w == 1: yield f, w, p, p + 8; p += 8
        elif w == 5: yield f, w, p, p + 4; p += 4
        else: raise ValueError("wire type %d" % w)
def decode_points(buf, z, r, c):
    b = memoryview(buf).tobytes(); span = 288 / 2 ** z; lon0 = -180 + c * span; lat0 = 90 - r * span; out = []
    for f, w, v, _ in _msg(b, 0, len(b)):
        if f != 3: continue
        keys, vals, feats, extent = [], [], [], 4096
        for f2, w2, v2, _ in _msg(b, v[0], v[1]):
            if f2 == 3: keys.append(b[v2[0]:v2[1]].decode())
            elif f2 == 5: extent = v2
            elif f2 == 4:
                val = None
                for f3, w3, v3, _ in _msg(b, v2[0], v2[1]):
                    if f3 == 1: val = b[v3[0]:v3[1]].decode()
                    elif f3 == 2: val = struct.unpack("<f", b[v3:v3 + 4])[0]
                    elif f3 == 3: val = struct.unpack("<d", b[v3:v3 + 8])[0]
                    elif f3 in (4, 5): val = v3
                    elif f3 == 6: val = (v3 >> 1) ^ -(v3 & 1)
                    elif f3 == 7: val = bool(v3)
                vals.append(val)
            elif f2 == 2: feats.append((v2[0], v2[1]))
        for a, e in feats:
            tags, geom = [], []
            for f3, w3, v3, _ in _msg(b, a, e):
                if w3 == 2 and f3 in (2, 4):
                    p, arr = v3[0], []
                    while p < v3[1]: x, p = _varint(b, p); arr.append(x)
                    if f3 == 2: tags = arr
                    else: geom = arr
            props = {keys[tags[i]]: vals[tags[i + 1]] for i in range(0, len(tags), 2)}
            if "LATITUDE" in props: out.append(props)
    return out
def fetch_day(sat, day):
    rows = []
    for r, c in TILES:
        url = GIBS.format(sat=sat, day=day, z=Z, r=r, c=c)
        try:
            with urllib.request.urlopen(url, timeout=60) as resp: buf = resp.read()
        except urllib.error.HTTPError as e:
            if e.code in (404, 204): continue
            raise
        if not buf: continue
        for p in decode_points(buf, Z, r, c):
            la, lo = float(p["LATITUDE"]), float(p["LONGITUDE"])
            if not (BBOX[0] <= la <= BBOX[1] and BBOX[2] <= lo <= BBOX[3]): continue
            rows.append([SAT_CODE[sat], str(p["ACQ_DATE"]).replace("-", ""), str(p["ACQ_TIME"]).replace(":", "").zfill(4),
                         round(la, 5), round(lo, 5), str(p["CONFIDENCE"])[0], float(p.get("FRP") or 0), p.get("DAYNIGHT", "D")])
    return rows

# ---------- assignment to Bangladesh, static sources, weights, land cover ----------
SEAS = {1: "NovFeb", 2: "NovFeb", 11: "NovFeb", 12: "NovFeb", 3: "Mar", 4: "Apr", 5: "May"}
def season(m): return SEAS.get(m, "JunOct")
class Processor:
    def __init__(self, ref_path, coef_path, lc_path):
        from scipy.spatial import cKDTree
        R = np.load(ref_path); bd = R["bd"].astype(bool); c = math.cos(math.radians(23.7))
        self.t_bd = cKDTree(np.column_stack([R["lat"][bd], R["lon"][bd] * c]))      # FIRMS Bangladesh detections 2012-2025
        self.t_fx = cKDTree(np.column_stack([R["lat"][~bd], R["lon"][~bd] * c]))    # neighbouring-country detections, March 2025
        self.static = set(map(tuple, R["static_cells"].tolist()))
        self.coef = json.load(open(coef_path))
        with gzip.open(lc_path, "rb") as f: self.lc = np.load(io.BytesIO(f.read()))
        self.lc_bbox = (20.4, 26.8, 87.8, 93.0)
    ALPHA = 4.0   # Bangladesh if 4x closer to past Bangladeshi fires than to foreign ones (held-out test: precision 0.985, recall 0.992)
    def in_bd(self, lat, lon):
        X = np.column_stack([lat, np.asarray(lon) * math.cos(math.radians(23.7))])
        db, _ = self.t_bd.query(X); df, _ = self.t_fx.query(X)
        return db * self.ALPHA < df, db
    def static_cell(self, lat, lon): return (int(round(lat / 0.004)), int(round(lon / 0.004)))
    def weight(self, src, lat, lon, month):
        if src == "VS": return 1.0
        key = "n20" if src == "V1" else "n21"; s = season(month)
        cell = int(math.floor(lat / 0.5)) * 1000 + int(math.floor(lon / 0.5))
        return float(self.coef[key]["cell"].get(f"{s}|{cell}", self.coef[key]["season"][s]))
    def landcover(self, lat, lon):
        S, N, W, E = self.lc_bbox; H, Wd = self.lc.shape
        r = min(H - 1, max(0, int((N - lat) / (N - S) * H))); c = min(Wd - 1, max(0, int((lon - W) / (E - W) * Wd)))
        return int(self.lc[r, c])
    def process(self, rows):
        if not rows: return []
        lat = np.array([r[3] for r in rows]); lon = np.array([r[4] for r in rows]); bd, dist = self.in_bd(lat, lon)
        out = []
        for r, ok, dd in zip(rows, bd, dist):
            if not ok: continue
            src, date, tm, la, lo, conf, frp, dn = r
            typ = 2 if self.static_cell(la, lo) in self.static else 0
            out.append([src, date, tm, la, lo, conf, round(frp, 2), dn, typ, round(self.weight(src, la, lo, int(date[4:6])), 3), self.landcover(la, lo)])
        return out
