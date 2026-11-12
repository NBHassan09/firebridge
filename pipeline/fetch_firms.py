"""Download NASA FIRMS country yearly summaries for Bangladesh and write data_in/firms_bd.csv.gz.
Source: https://firms.modaps.eosdis.nasa.gov/country/  (standard processing; no key needed)
Usage: python pipeline/fetch_firms.py"""
import csv, gzip, io, os, urllib.request
BASE = "https://firms.modaps.eosdis.nasa.gov/data/country/{s}/{y}/{s}_{y}_Bangladesh.csv"
JOBS = [("modis", y) for y in range(2000, 2026)] + [("viirs-snpp", y) for y in range(2012, 2026)] + [("viirs-jpss1", y) for y in range(2018, 2026)]
os.makedirs("data_in", exist_ok=True)
out = io.StringIO(); w = csv.writer(out); w.writerow(["src","date","time","lat","lon","conf","frp","type","daynight","version"])
for s, y in JOBS:
    try:
        txt = urllib.request.urlopen(BASE.format(s=s, y=y), timeout=120).read().decode()
    except Exception as e:
        print("skip", s, y, e); continue
    for r in csv.DictReader(io.StringIO(txt)):
        src = ("MT" if r["satellite"].startswith("T") else "MA") if s == "modis" else ("VS" if s == "viirs-snpp" else "V1")
        w.writerow([src, r["acq_date"].replace("-", ""), r["acq_time"].zfill(4), f'{float(r["latitude"]):.4f}', f'{float(r["longitude"]):.4f}',
                    r["confidence"], r["frp"], r["type"], r["daynight"], r["version"]])
    print("ok", s, y)
with gzip.open("data_in/firms_bd.csv.gz", "wt") as f: f.write(out.getvalue())
