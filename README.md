# FireBridge

A harmonized burning-activity calendar for Bangladesh, built for the 2026 NASA Space Apps Challenge "Harmonization of MODIS and VIIRS Hot Spots".

Pick an area (draw a box or shape, enter a point and radius, or choose an example) and a date range. The app shows:
- a months × years calendar of hotspot counts from 2001 to 2025, coloured by count or by rank against the same month in other years;
- how a selected month compares with its historical median, with an 80% uncertainty range for months estimated from MODIS;
- which land cover burned each year, with a land-cover filter for the whole page;
- the weather and air in the selected month against normal, and a test of whether rainfall, humidity, temperature or soil wetness explain the fire count across years;
- a check of the harmonization inside the chosen area (MODIS-based estimate vs VIIRS S-NPP observed, plus the NOAA-20 link);
- CSV downloads of the hotspots and the monthly calendar.

## Data (all NASA)
- Hotspots: NASA FIRMS country yearly summaries, standard processing. MODIS C6.1 (MCD14ML, Terra and Aqua), VIIRS S-NPP 375 m (VNP14IMG), VIIRS NOAA-20 375 m (VJ114IMG). Nov 2000 to Dec 2025, 279,878 detections.
- Land cover: MODIS MCD12Q1 IGBP land cover, annual 2001–2024, via NASA GIBS. Each hotspot takes the class of its pixel in its own year.
- Weather and air: NASA POWER monthly regional API, 2001–2025: precipitation (PRECTOTCORR), relative humidity (RH2M), maximum temperature (T2M_MAX), topsoil wetness (GWETTOP), aerosol optical depth (AOD_55).
- Basemap: NASA GIBS WMS (Blue Marble shaded relief, VIIRS Black Marble). Border lines: US Department of State LSIB, served by NASA GIBS.

## Live data, no keys
- 2026 onward: VIIRS S-NPP, NOAA-20 and NOAA-21 detections from NASA GIBS vector tiles (`VIIRS_*_Thermal_Anomalies_375m_All`, the FIRMS near-real-time feed), fetched daily by a GitHub Action (`.github/workflows/update-live.yml`, `pipeline/update_live.py`). NASA POWER weather for the current year is fetched the same way.
- Assignment to Bangladesh: a detection counts when it is four times closer to past FIRMS Bangladesh detections than to neighbouring-country ones. Held-out March 2025 days: recall 99.2%, precision 98.5%, monthly count within 1% of the FIRMS country file. Every one of the 7,616 FIRMS Bangladesh detections in March 2025 was present in the GIBS feed.
- Handover: S-NPP ended on 2 Nov 2026. From November 2026 the series is NOAA-20 converted to S-NPP equivalents. Jan–Oct 2026 check: converted NOAA-20 = 94% of S-NPP (inside the 2019–2025 range of 94–109%); March 2026 was the weakest month at 84%.
- A month in progress is compared with the same days of that month in past years.

## Method
1. Filter: type 0 (presumed vegetation fire); MODIS confidence ≥ 30; VIIRS nominal or high. Static land sources (type 2: kilns, flares) and offshore (type 3) are excluded from counts but exported.
2. Reference: VIIRS S-NPP from February 2012.
3. Convert MODIS: each MODIS detection gets weight k = S-NPP/MODIS ratio for its season (Nov–Feb, Mar, Apr, May, Jun–Oct) and 0.5° cell over 2012–2021, shrunk toward the season ratio (α = 20) and rescaled to keep season totals. Terra-only months (to Jun 2002) use a Terra-only ratio and are flagged low confidence.
4. Continuity: S-NPP delivery ended 2 Nov 2026; NOAA-20 is linked to S-NPP with the same estimator (2019–2025).
5. Uncertainty: 80% ranges from the 10th/90th percentiles of observed/estimated ratios in leave-one-year-out testing, binned by count (83% coverage in testing).

## What the context adds
- What burned: 52% of VIIRS hotspots fall on shrub, savanna and grass (mostly woody savanna, i.e. jhum fallow in the Chittagong Hill Tracts), 31% on forest, 15% on cropland.
- Why some years burn more: across 2003–2025, March–April hotspots nationally fall as January–April rainfall rises (Spearman ρ = −0.71) and as March–April humidity rises (ρ = −0.66); the links hold within VIIRS-only years (ρ = −0.62 and −0.73), so they are not an artefact of the MODIS conversion. The app runs the same test, with a permutation p-value, for any drawn area and month.
- Aerosol optical depth over Bangladesh rose while fires fell, so it is shown as air-quality context, not as a smoke signal.

## Early warning
- Season tracker: cumulative hotspots 1 Jan–31 May against the P10–P90 range of past years at the same date.
- What the data supports: rain shifts burning in time. After a wetter January–March, about 67% of March–April burning came in April, against 46% after a drier start (ρ = −0.46 nationally, −0.53 in the Hill Tracts). Early-season counts do not predict the season total (ρ ≈ 0), and pre-season (Nov–Feb) rain is only a weak predictor of the total (tercile hit rate 52% vs 33% by chance), so the app reports timing, not a total forecast.

## Limitations
- Small, smouldering fires under the Sundarbans mangrove canopy are mostly missed by 375 m and 1 km sensors (e.g. no detections during the reported March 2025 fires).
- Terra-only months (Nov 2000–Jun 2002) are low confidence; 2022 MODIS Aqua under-recorded spring hotspots.
- Near-real-time 2026 data will be replaced by FIRMS standard processing a few months later.
- MERRA-2 aerosol over Bangladesh is dominated by regional pollution, so it is context, not a smoke measure.

## Validation (leave-one-year-out, 2012–2021)
| Method | National monthly error | Cell-month error | Cell annual rank ρ |
|---|---|---|---|
| One national ratio | 44% | 38% | 0.91 |
| Ratio by season | 20% | 37% | 0.92 |
| Season × 0.5° cell, shrunk (used) | 18% | 36% | 0.93 |

Out of sample, the MODIS-based estimate is within 3–7% of S-NPP in 2023–2025 despite orbital drift of both MODIS platforms. In 2022 MODIS Aqua recorded about 60% of the expected spring hotspots. The NOAA-20 to S-NPP link has an 8% national monthly error.

### Control test: was the 2012 jump real?

Aqua and Terra MODIS flew unchanged across the 2012 boundary, so a MODIS-only series
is an instrument that did not change. If burning had really surged, MODIS alone would
show it. Mean hotspots per year, 2003–2011 against 2013–2021 (the 2012 boundary year
is excluded):

| Series | 2003–2011 | 2013–2021 | Change |
|---|---|---|---|
| Aqua MODIS only (control) | 2,922 | 2,041 | **−30.1%** |
| Terra MODIS only (control) | 285 | 234 | −17.8% |
| MODIS Terra + Aqua | 3,207 | 2,276 | −29.0% |
| Raw splice: MODIS, then VIIRS | 3,207 | 9,157 | **+185.5%** |
| Harmonized (FireBridge) | 12,514 | 9,157 | **−26.8%** |

Burning in Bangladesh **fell** by about 30% between those periods. The raw splice
reports a 186% rise — wrong by 216 percentage points, and in the opposite direction.
The harmonized series lands within **3.3 percentage points** of the unchanged-instrument
control. Aqua is the primary control: its ~13:30 overpass is closest to S-NPP and its
orbit was stable through 2021. Terra is a second control; its overpass time began
drifting after 2020, which is why the comparison stops there.

Run it with `python pipeline/control.py` (writes `pipeline/control_test.json`).

## Reproduce
```
python pipeline/fetch_firms.py        # downloads FIRMS CSVs -> data_in/firms_bd.csv.gz
mkdir -p analysis app && cp pipeline/{common,harmonize,export}.py analysis/
python analysis/harmonize.py          # fits, validates, writes analysis/coefficients.json
python analysis/context.py            # needs data_in/nasa_context.json (POWER + land cover); writes app/context.json
python analysis/export.py             # writes app/meta.json and app/points.bin
python pipeline/update_live.py        # 2026+ near-real-time detections and weather (no keys)
python -c "import gzip,base64;open('app/points.gz.b64.txt','wb').write(base64.b64encode(gzip.compress(open('app/points.bin','rb').read(),9)))"
```
The site is static: serve the files in `site/` (Vercel: import the folder, no build step).
