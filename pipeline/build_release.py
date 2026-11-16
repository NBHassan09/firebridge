"""Build the public data release for Zenodo.

Produces a self-contained, documented archive of the harmonized Bangladesh
active-fire record: the per-detection table (with its harmonization weight),
the national monthly series with 80% uncertainty intervals, the 0.5-degree
cell-month series, the fitted coefficients, and the two validation experiments.

Usage:  python pipeline/build_release.py            (run from the project root)
Writes: release/ and release/firebridge-bangladesh-v1.0.zip
"""
import json, os, shutil, zipfile
import numpy as np
import pandas as pd

VERSION = '1.0'
OUT = 'release'
IGBP = {0: 'water', 1: 'evergreen needleleaf forest', 2: 'evergreen broadleaf forest',
        3: 'deciduous needleleaf forest', 4: 'deciduous broadleaf forest', 5: 'mixed forest',
        6: 'closed shrubland', 7: 'open shrubland', 8: 'woody savanna', 9: 'savanna',
        10: 'grassland', 11: 'permanent wetland', 12: 'cropland', 13: 'urban and built-up',
        14: 'cropland/natural mosaic', 15: 'snow and ice', 16: 'barren', 255: 'unclassified'}
SENSOR = {'MT': 'MODIS Terra', 'MA': 'MODIS Aqua',
          'VS': 'VIIRS S-NPP', 'V1': 'VIIRS NOAA-20'}
ERA = {0: 'terra-only (low confidence)', 1: 'MODIS-harmonized', 2: 'VIIRS S-NPP observed'}


def interval(est, bins, q):
    """80% interval for an estimated count, from the empirical ratio percentiles."""
    b = np.clip(np.digitize(est, bins[1:], right=False), 0, len(q) - 1)
    lo = est * np.array([q[i][0] for i in b])
    hi = est * np.array([q[i][1] for i in b])
    return lo.round(1), hi.round(1)


def main():
    os.makedirs(OUT, exist_ok=True)
    d = pd.read_pickle('analysis/detections.pkl')
    d = d.reset_index(drop=True)
    d['igbp'] = pd.read_pickle('analysis/igbp.pkl').reset_index(drop=True).igbp.values
    coef = json.load(open('analysis/coefficients.json'))

    # ---------- 1. per-detection table ----------
    # `passes_filter` = a vegetation fire at usable confidence.
    # `in_series`     = actually contributes to the harmonized series.  These differ:
    #   NOAA-20 detections pass the filter but carry weight 0, because S-NPP is the
    #   reference sensor and counting both during their overlap would double-count.
    passes = (d.confok & (d.type == 0)).values
    det = pd.DataFrame({
        'date': d.dt.dt.strftime('%Y-%m-%d'),
        'time_utc': d.time.astype(str).str.zfill(4),
        'latitude': d.lat.round(5),
        'longitude': d.lon.round(5),
        'sensor': d.src.map(SENSOR),
        'confidence': d.conf,
        'frp_mw': d.frp,
        'fire_type': d.type,
        'daynight': d.daynight,
        'passes_filter': passes,
        'in_series': passes & (d.w.values > 0),
        'weight': d.w.round(4),
        'era': d.era.map(ERA),
        'landcover_igbp': d.igbp,
        'landcover_name': d.igbp.map(IGBP).fillna('unclassified'),
    }).sort_values(['date', 'time_utc']).reset_index(drop=True)
    det.to_csv(f'{OUT}/firebridge_detections_v{VERSION}.csv.gz', index=False,
               compression={'method': 'gzip', 'mtime': 0})

    # ---------- 2. national monthly series ----------
    f = d[d.confok & (d.type == 0)]
    U = coef['uncertainty']['modis_era']
    Ut = coef['uncertainty']['terra_only_era']
    nat = f.groupby('ym').agg(harmonized=('w', 'sum'), era=('era', 'max')).reset_index()
    nat['year'] = nat.ym // 100
    nat['month'] = nat.ym % 100
    nat['harmonized'] = nat.harmonized.round(1)
    obs = f[f.src == 'VS'].groupby('ym').size()
    nat['snpp_observed'] = nat.ym.map(obs)
    lo_m, hi_m = interval(nat.harmonized.values, U['bins'] + [1e9], U['q10_q90_ratio'])
    lo_t, hi_t = interval(nat.harmonized.values, Ut['bins'] + [1e9], Ut['q10_q90_ratio'])
    terra = (nat.era == 0).values
    viirs = (nat.era == 2).values
    nat['lower_80'] = np.where(terra, lo_t, lo_m)
    nat['upper_80'] = np.where(terra, hi_t, hi_m)
    # observed years carry no estimation interval
    nat.loc[viirs, ['lower_80', 'upper_80']] = np.nan
    # The Terra-only era (to Jun 2002) needs weights of 38-60 applied to a few hundred
    # detections.  It is published for completeness but must not be used for analysis.
    nat['low_confidence'] = terra
    nat['era'] = nat.era.map(ERA)
    nat[['year', 'month', 'harmonized', 'lower_80', 'upper_80',
         'snpp_observed', 'era', 'low_confidence']].to_csv(
             f'{OUT}/firebridge_monthly_national_v{VERSION}.csv', index=False)

    # ---------- 3. cell-month series ----------
    g = f.copy()
    g['cell_lat'] = np.floor(g.lat / 0.5) * 0.5
    g['cell_lon'] = np.floor(g.lon / 0.5) * 0.5
    cm = (g.groupby(['y', 'm', 'cell_lat', 'cell_lon'])
            .agg(harmonized=('w', 'sum'), era=('era', 'max')).reset_index())
    cm = cm[cm.harmonized > 0]
    lo_m, hi_m = interval(cm.harmonized.values, U['bins'] + [1e9], U['q10_q90_ratio'])
    cm['lower_80'] = lo_m
    cm['upper_80'] = hi_m
    cm.loc[cm.era == 2, ['lower_80', 'upper_80']] = np.nan
    cm['harmonized'] = cm.harmonized.round(1)
    cm['low_confidence'] = (cm.era == 0).values
    cm['era'] = cm.era.map(ERA)
    cm = cm.rename(columns={'y': 'year', 'm': 'month'})
    cm[['year', 'month', 'cell_lat', 'cell_lon', 'harmonized',
        'lower_80', 'upper_80', 'era', 'low_confidence']].to_csv(
            f'{OUT}/firebridge_monthly_cells_v{VERSION}.csv.gz', index=False,
            compression={'method': 'gzip', 'mtime': 0})

    # ---------- 4. coefficients and validation ----------
    for src, dst in [('analysis/coefficients.json', 'firebridge_coefficients'),
                     ('analysis/control_test.json', 'firebridge_control_test'),
                     ('analysis/travel_test.json', 'firebridge_transferability_test')]:
        shutil.copy(src, f'{OUT}/{dst}_v{VERSION}.json')

    # ---------- 5. manifest ----------
    rows = []
    for fn in sorted(os.listdir(OUT)):
        if fn.endswith('.zip'):
            continue
        rows.append({'file': fn, 'bytes': os.path.getsize(f'{OUT}/{fn}')})
    json.dump({'dataset': 'FireBridge harmonized active-fire record for Bangladesh',
               'version': VERSION,
               'licence': 'CC BY 4.0',
               'archive_coverage': '2000-11-02 to 2025-12-31 (FIRMS standard processing)',
               'series_coverage': 'monthly, 2001-01 to 2025-12',
               'recommended_start': '2002-07 (earlier months are Terra-only, '
                                    'flagged low_confidence, and should not be used '
                                    'for trend analysis)',
               'detections_in_archive': int(len(det)),
               'detections_passing_filter': int(det.passes_filter.sum()),
               'detections_in_series': int(det.in_series.sum()),
               'files': rows},
              open(f'{OUT}/manifest.json', 'w'), indent=1)

    print(f'rows: detections {len(det):,}  national months {len(nat):,}  cell-months {len(cm):,}')
    print(f'passes filter {det.passes_filter.sum():,}   in series {det.in_series.sum():,}')
    for r in rows:
        print(f"  {r['file']:52s} {r['bytes']/1e6:8.2f} MB")


if __name__ == '__main__':
    main()
