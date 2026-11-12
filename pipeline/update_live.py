"""Daily update for FireBridge, no API keys needed.
- VIIRS NOAA-20 and NOAA-21 fire detections from NASA GIBS vector tiles (S-NPP ended 2 Nov 2026)
- NASA POWER daily weather, aggregated to monthly means
Writes site/live.json and site/live_power.json. Run from the repo root: python pipeline/update_live.py"""
import json, os, sys, datetime as dt, urllib.request
sys.path.insert(0, os.path.dirname(__file__))
from live_common import fetch_day, Processor

SITE = "site"; REF = "pipeline/live_ref.npz"; COEF = "pipeline/link_coeffs.json"; LC = "pipeline/landcover_2024.npy.gz"
def load(path, default):
    try: return json.load(open(path))
    except FileNotFoundError: return default
live = load(f"{SITE}/live.json", {"fields": [], "rows": []})
proc = Processor(REF, COEF, LC)
today = dt.date.today(); start = today - dt.timedelta(days=4)               # re-fetch a few days: NRT files fill in late
keep = [r for r in live["rows"] if r[1] < start.strftime("%Y%m%d")]
new = []; failed = False
for k in range((today - start).days):
    day = (start + dt.timedelta(days=k)).isoformat()
    for sat in ("NOAA20", "NOAA21"):
        try: new += proc.process(fetch_day(sat, day))
        except Exception as e: print("skip", sat, day, e); failed = True
if failed:   # never drop days we already had because a fetch failed
    have = {r[1] for r in new}
    keep += [r for r in live["rows"] if r[1] >= start.strftime("%Y%m%d") and r[1] not in have]
seen = set(); rows = []
for r in keep + new:
    key = (r[0], r[1], r[2], r[3], r[4])
    if key in seen: continue
    seen.add(key); rows.append(r)
rows.sort(key=lambda r: (r[1], r[2]))
live.update(rows=rows, updated=today.isoformat(), last_day=rows[-1][1] if rows else live.get("last_day"))
json.dump(live, open(f"{SITE}/live.json", "w"), separators=(",", ":"))
print("live rows", len(rows), "new", len(new))

# NASA POWER, current year, daily -> monthly means
power = load(f"{SITE}/live_power.json", {"params": {}})
y = today.year
for p in ("PRECTOTCORR", "RH2M", "T2M_MAX", "GWETTOP", "AOD_55"):
    url = ("https://power.larc.nasa.gov/api/temporal/daily/regional?parameters=%s&community=AG&latitude-min=20.5&latitude-max=26.75"
           "&longitude-min=88&longitude-max=92.75&start=%d0101&end=%s&format=JSON") % (p, y, (today - dt.timedelta(days=2)).strftime("%Y%m%d"))
    try:
        j = json.load(urllib.request.urlopen(url, timeout=120)); fill = j.get("header", {}).get("fill_value", -999); pts = []
        for f in j["features"]:
            agg = {}
            for k, x in f["properties"]["parameter"][p].items():
                if x is None or x == fill or x <= -998: continue
                a = agg.setdefault(k[:6], [0.0, 0]); a[0] += x; a[1] += 1
            pts.append({"lon": f["geometry"]["coordinates"][0], "lat": f["geometry"]["coordinates"][1],
                        "m": {m: [round(s / n, 3), n] for m, (s, n) in agg.items()}})
        power["params"][p] = pts
    except Exception as e: print("power skip", p, e)
power["source"] = "NASA POWER daily regional API, aggregated to monthly means; updated %s" % today.isoformat()
json.dump(power, open(f"{SITE}/live_power.json", "w"), separators=(",", ":"))
