# FireBridge

A harmonized burning-activity calendar for Bangladesh, built for the 2026 NASA Space Apps Challenge "Harmonization of MODIS and VIIRS Hot Spots".

Pick an area (draw a box or shape, enter a point and radius, or choose an example) and a date range. The app shows:
- a months × years calendar of hotspot counts from 2001 to 2025, coloured by count or by rank against the same month in other years;
- how a selected month compares with its historical median, with an 80% uncertainty range for months estimated from MODIS;
- a check of the harmonization inside the chosen area (MODIS-based estimate vs VIIRS S-NPP observed, plus the NOAA-20 link);
- CSV downloads of the hotspots and the monthly calendar.

## Data (all NASA)
- Hotspots: NASA FIRMS country yearly summaries, standard processing. MODIS C6.1 (MCD14ML, Terra and Aqua), VIIRS S-NPP 375 m (VNP14IMG), VIIRS NOAA-20 375 m (VJ114IMG). Nov 2000 to Dec 2025, 279,878 detections.
- Basemap: NASA GIBS WMS (Blue Marble shaded relief, VIIRS Black Marble). Border lines: US Department of State LSIB, served by NASA GIBS.

## Method
1. Filter: type 0 (presumed vegetation fire); MODIS confidence ≥ 30; VIIRS nominal or high. Static land sources (type 2: kilns, flares) and offshore (type 3) are excluded from counts but exported.
2. Reference: VIIRS S-NPP from February 2012.
3. Convert MODIS: each MODIS detection gets weight k = S-NPP/MODIS ratio for its season (Nov–Feb, Mar, Apr, May, Jun–Oct) and 0.5° cell over 2012–2021, shrunk toward the season ratio (α = 20) and rescaled to keep season totals. Terra-only months (to Jun 2002) use a Terra-only ratio and are flagged low confidence.
4. Continuity: S-NPP delivery ended 2 Nov 2026; NOAA-20 is linked to S-NPP with the same estimator (2019–2025).
5. Uncertainty: 80% ranges from the 10th/90th percentiles of observed/estimated ratios in leave-one-year-out testing, binned by count (83% coverage in testing).

## Validation (leave-one-year-out, 2012–2021)
| Method | National monthly error | Cell-month error | Cell annual rank ρ |
|---|---|---|---|
| One national ratio | 44% | 38% | 0.91 |
| Ratio by season | 20% | 37% | 0.92 |
| Season × 0.5° cell, shrunk (used) | 18% | 36% | 0.93 |

Out of sample, the MODIS-based estimate is within 3–7% of S-NPP in 2023–2025 despite orbital drift of both MODIS platforms. In 2022 MODIS Aqua recorded about 60% of the expected spring hotspots. The NOAA-20 to S-NPP link has an 8% national monthly error.

## Reproduce
```
python pipeline/fetch_firms.py        # downloads FIRMS CSVs -> data_in/firms_bd.csv.gz
mkdir -p analysis app && cp pipeline/{common,harmonize,export}.py analysis/
python analysis/harmonize.py          # fits, validates, writes analysis/coefficients.json
python analysis/export.py             # writes app/meta.json and app/points.bin
python -c "import gzip,base64;open('app/points.gz.b64.txt','wb').write(base64.b64encode(gzip.compress(open('app/points.bin','rb').read(),9)))"
```
The site is static: serve the files in `site/` (Vercel: import the folder, no build step).
