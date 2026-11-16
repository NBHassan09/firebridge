# FireBridge: a harmonized MODIS–VIIRS active fire record for Bangladesh, v1.0

A single, internally comparable monthly record of vegetation burning in Bangladesh
from January 2001 to December 2025, built by converting NASA FIRMS MODIS detections
into VIIRS S-NPP equivalents.

**Why this exists.** MODIS (1 km) has observed Bangladesh since 2000 and VIIRS (375 m)
since February 2012, and VIIRS records about 3.9 vegetation-fire detections for every
one MODIS records. A time series assembled by using MODIS before 2012 and VIIRS after
it therefore measures the change of instrument, not the change in burning. Over
2003–2011 against 2013–2021 such a splice reports national burning rising 186%, while
an unchanged instrument (Aqua MODIS alone) records it falling 30%. This dataset
removes that discontinuity.

- **Licence:** CC BY 4.0
- **Source data:** NASA FIRMS country archives, standard processing
- **Code:** https://github.com/NBHassan09/firebridge
- **Interactive viewer:** https://firebridge-bangladesh.vercel.app

---

## Files

| File | Rows | What it is |
|---|---|---|
| `firebridge_detections_v1.0.csv.gz` | 279,878 | Every detection in the archive, with its harmonization weight |
| `firebridge_monthly_national_v1.0.csv` | 278 | National monthly harmonized counts with 80% intervals |
| `firebridge_monthly_cells_v1.0.csv.gz` | 5,347 | The same, per 0.5° cell and month |
| `firebridge_coefficients_v1.0.json` | — | Fitted weights, season ratios, validation metrics |
| `firebridge_control_test_v1.0.json` | — | The unchanged-instrument control experiment |
| `firebridge_transferability_test_v1.0.json` | — | The same method applied to Portugal and Nepal |
| `manifest.json` | — | Coverage, counts, file sizes |

## Start here

Most users want `firebridge_monthly_national_v1.0.csv`:

```
year,month,harmonized,lower_80,upper_80,snpp_observed,era,low_confidence
2003,3,7735.6,5494.3,9134.1,,MODIS-harmonized,False
2012,3,5455.0,,,5455.0,VIIRS S-NPP observed,False
```

- `harmonized` — the VIIRS-equivalent count. **This is the column to use.**
- `lower_80` / `upper_80` — empirical 80% interval. Blank where the month is observed
  directly by VIIRS and so carries no estimation uncertainty.
- `snpp_observed` — raw S-NPP count, where available. Populated from Feb 2012.
- `era` — how this month was produced (see below).
- `low_confidence` — **drop these rows.** See the warning below.

## Three eras

| Era | Period | How the count is produced |
|---|---|---|
| `terra-only (low confidence)` | Jan 2001 – Jun 2002 | Terra alone, weights of 38–60. **Do not use.** |
| `MODIS-harmonized` | Jul 2002 – Jan 2012 | Terra + Aqua converted to S-NPP equivalents |
| `VIIRS S-NPP observed` | Feb 2012 – Dec 2025 | S-NPP itself, weight 1. Not modelled. |

> **Warning on the Terra-only era.** Before Aqua's first Bangladesh detection in
> September 2002, only Terra observed the country, and its mid-morning overpass misses
> most burning. Reconstructing S-NPP-equivalent counts from it needs weights of 38 to
> 60 applied to a few hundred detections, which cross-validates at 50% national monthly
> error against 18% for the two-platform model. These 11 months are published so the
> archive is complete, not because they are usable. **For any trend analysis, start at
> July 2002, and preferably at 2003.**

## Per-detection table

Two boolean columns are easy to confuse, and using the wrong one will double-count:

- `passes_filter` (247,320 rows) — a presumed vegetation fire at usable confidence.
- `in_series` (152,351 rows) — **actually contributes to the harmonized series.**

They differ because NOAA-20 detections pass the filter but carry weight 0: S-NPP is the
reference sensor, and counting both during their 2018–2025 overlap would double-count
the same fires. NOAA-20 is retained in the table because it is what continues the record
after S-NPP ceased delivery on 2 November 2026.

**To reproduce the national series:** sum `weight` over rows where `in_series` is true.
Do not count rows.

Columns: `date`, `time_utc` (HHMM, zero-padded), `latitude`, `longitude`, `sensor`,
`confidence`, `frp_mw`, `fire_type`, `daynight`, `passes_filter`, `in_series`, `weight`,
`era`, `landcover_igbp`, `landcover_name`.

## Filtering applied

- **Fire type:** only type 0 (presumed vegetation fire) is counted. Type 2 (static land
  sources — in Bangladesh overwhelmingly brick kilns and gas flares, 6,054 records) and
  type 3 (offshore, 1,336 records) are retained in the detection table but excluded from
  all counts. Brick kilns are persistent, bright and clustered around Dhaka; including
  them would both inflate counts and invent a spatial pattern.
- **Confidence:** MODIS at ≥ 30; VIIRS at nominal or high. The two scales are not
  interchangeable, so a single numeric threshold is not applied to both.

## Method in brief

Each MODIS detection gets a weight

    k(s,c) = (Y(s,c) + α·k̂(s)) / (X(s,c) + α)

where X and Y are MODIS and S-NPP detection counts in season `s` and 0.5° cell `c` over
the 2012–2021 overlap, k̂(s) is the season-wide ratio, and α = 20 MODIS detections.
Weights are then rescaled per season to preserve observed season totals.

Seasons are Nov–Feb, Mar, Apr, May, Jun–Oct. March and April are separate single-month
seasons because they carry 96% of national burning. Season ratios: 8.23 (Nov–Feb),
3.48 (Mar), 3.64 (Apr), 5.10 (May), 3.46 (Jun–Oct).

## Accuracy

Leave-one-year-out cross-validation over 2012–2021, against held-out S-NPP observations:

| Estimator | National monthly error | Cell-month WAPE | Cell annual rank ρ |
|---|---|---|---|
| Single national ratio | 43.7% | 37.9% | 0.907 |
| Ratio by season | 20.1% | 37.2% | 0.915 |
| **Season × cell, shrunk (this dataset)** | **18.0%** | **35.5%** | **0.928** |

Applied to 2022–2025, years entirely outside the training window, the MODIS-based
estimate lands within 3–7% of observed S-NPP in 2023, 2024 and 2025. In 2022 it reaches
only 60%, because MODIS Aqua itself under-recorded that spring — a failure of the input
data, not of the conversion.

Observed values fall inside the published 80% intervals 83.2% of the time.

## Known limitations

1. **Small fires under closed canopy are missed by both sensors.** No detections were
   recorded during the reported March 2025 Sundarbans mangrove fires. Absence of
   detections is not absence of fire.
2. **2022 MODIS Aqua under-recorded** spring detections by roughly 40%.
3. **Orbital drift.** Terra's median local overpass over Bangladesh moved from 10.70 in
   2021 to 9.70 by 2025, Aqua's from 13.33 to 14.20. Weights were fitted on a stable
   period; expect MODIS-based estimates to degrade as drift continues.
4. **Single cell, single month, MODIS era** is the weakest cell of the product. The
   intervals are wide there for a reason.
5. **The spatial refinement is regime-specific.** Splitting the ratio by cell helps in
   Bangladesh, where jhum fires recur in the same fields, and hurts in Portugal, where
   wildfires move. See the transferability file before reusing the method elsewhere.

## Citation

> Hassan, N. B. (2026). FireBridge: a harmonized MODIS–VIIRS active fire record for
> Bangladesh, 2001–2025 (v1.0) [Data set]. Zenodo.

## Acknowledgement of source data

Produced using NASA FIRMS data (MODIS Collection 6.1 MCD14ML; VIIRS VNP14IMG and
VJ114IMG), MODIS MCD12Q1 v6.1 land cover, and NASA POWER meteorology. FIRMS is operated
by NASA's Earth Science Data and Information System. The author is not affiliated with
NASA, and this dataset is not a NASA product.
