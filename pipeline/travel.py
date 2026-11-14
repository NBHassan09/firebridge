"""Does the FireBridge harmonization travel to other countries?

The Bangladesh model groups months into five hand-chosen seasons.  That grouping is
Bangladesh-specific, so for a fair cross-country test every country here is run
through one identical, country-agnostic configuration:

    weight k(calendar month, 0.5-degree cell), empirical-Bayes shrinkage alpha=20
    toward the month-wide ratio, fitted on the 2012-2021 MODIS/S-NPP overlap,
    leave-one-year-out cross-validation.

Nothing is tuned per country.  Bangladesh is re-run under the same configuration so
the comparison is like for like.

Usage:  python pipeline/travel.py <extra_countries.csv.gz>
"""
import sys, json
import pandas as pd, numpy as np
from scipy.stats import spearmanr

ALPHA = 20; RES = 0.5
TRAIN = list(range(2012, 2022))
SWITCH = 201202              # first full month of the S-NPP record

def prep(d):
    d = d.copy()
    d['dt'] = pd.to_datetime(d.date.astype(str), format='%Y%m%d')
    d['y'] = d.dt.dt.year; d['m'] = d.dt.dt.month; d['ym'] = d.y * 100 + d.m
    isM = d.src.isin(['MT', 'MA'])
    d['confok'] = np.where(isM, pd.to_numeric(d.conf, errors='coerce') >= 30, d.conf.astype(str) != 'l')
    d['s'] = d.m                                    # season = calendar month (country-agnostic)
    d['cell'] = np.floor(d.lat / RES).astype(int) * 1000 + np.floor(d.lon / RES).astype(int)
    return d[(d.type == 0) & d.confok]

def fit_k(train, pred, ref='VS'):
    X = train[train.src.isin(pred)].groupby(['s', 'cell']).size().rename('X')
    Y = train[train.src == ref].groupby(['s', 'cell']).size().rename('Y')
    J = pd.concat([X, Y], axis=1).fillna(0)
    ks = J.groupby(level=0).Y.sum() / J.groupby(level=0).X.sum()
    J['k0'] = ks.reindex(J.index.get_level_values(0)).values
    J['k'] = (J.Y + ALPHA * J.k0) / (J.X + ALPHA)
    fac = J.groupby(level=0).Y.sum() / (J.k * J.X).groupby(level=0).sum()
    J['k'] = J.k * fac.reindex(J.index.get_level_values(0)).values
    return J, ks

def apply_k(df, J, ks):
    kd = J.k.to_dict(); kf = ks.to_dict(); g = float(np.nanmedian(list(kf.values())))
    return np.array([kd.get((s, c), kf.get(s, g)) for s, c in zip(df.s, df.cell)])

def loyo(ov, pred, mode='cell'):
    rows = []
    for yt in TRAIN:
        tr, te = ov[ov.y != yt], ov[ov.y == yt]
        if len(te) == 0 or len(tr[tr.src == 'VS']) == 0:
            continue
        xm = te[te.src.isin(pred)].copy()
        if len(xm) == 0:
            continue
        if mode == 'cell':
            J, ks = fit_k(tr, pred); xm['w'] = apply_k(xm, J, ks)
        elif mode == 'month':
            ks = tr[tr.src == 'VS'].groupby('s').size() / tr[tr.src.isin(pred)].groupby('s').size()
            g = float(ks.median()); xm['w'] = xm.s.map(ks).fillna(g).values
        else:                                       # one ratio for the whole country
            k = len(tr[tr.src == 'VS']) / max(len(tr[tr.src.isin(pred)]), 1)
            xm['w'] = k
        P = xm.groupby(['m', 'cell']).w.sum(); O = te[te.src == 'VS'].groupby(['m', 'cell']).size()
        R = pd.concat([P.rename('P'), O.rename('O')], axis=1).fillna(0); R['y'] = yt
        rows.append(R)
    return pd.concat(rows) if rows else None

def metrics(R):
    if R is None: return None
    nat = R.groupby(['y', 'm']).sum(); nb = nat[nat.O >= 50]
    ann = R.groupby(['y', 'cell']).sum(); ann = ann[ann.O >= 30]
    out = {'national_month_median_abs_pct_err': round(float((abs(nb.P - nb.O) / nb.O).median() * 100), 1),
           'cell_month_WAPE_pct': round(float(abs(R.P - R.O).sum() / R.O.sum() * 100), 1),
           'total_bias_pct': round(float((R.P.sum() / R.O.sum() - 1) * 100), 1),
           'months_tested': int(len(nb))}
    out['cell_annual_spearman'] = round(float(spearmanr(ann.P, ann.O).statistic), 3) if len(ann) > 10 else None
    return out

def run(name, f):
    ov = f[f.y.isin(TRAIN) & (f.ym >= SWITCH)]
    nV = int((ov.src == 'VS').sum()); nM = int(ov.src.isin(['MT', 'MA']).sum())
    yv = set(ov[ov.src == 'VS'].y.unique()); ym = set(ov[ov.src.isin(['MT', 'MA'])].y.unique())
    both = [int(x) for x in sorted(yv & ym)]
    if nV < 500 or nM < 200:
        return {'country': name, 'skipped': f'too few overlap detections (VIIRS {nV}, MODIS {nM})'}
    if len(both) < 6:
        return {'country': name, 'skipped': f'only {len(both)} complete overlap years retrieved '
                f'({both}); NASA FIRMS timed out on the rest, so no cross-validated figure is reported'}
    return {'country': name, 'complete_overlap_years': both,
            'overlap_detections': {'VIIRS_SNPP': nV, 'MODIS': nM, 'raw_VIIRS_per_MODIS': round(nV / nM, 2)},
            'cells_05deg': int(ov.cell.nunique()),
            'one_ratio': metrics(loyo(ov, ['MA', 'MT'], 'whole')),
            'month_ratio': metrics(loyo(ov, ['MA', 'MT'], 'month')),
            'month_x_cell_shrunk': metrics(loyo(ov, ['MA', 'MT'], 'cell'))}

if __name__ == '__main__':
    res = []
    bd = pd.read_csv('data_in/firms_bd.csv.gz', dtype={'time': str, 'conf': str, 'date': str}, low_memory=False)
    res.append(run('Bangladesh', prep(bd)))
    if len(sys.argv) > 1:
        ex = pd.read_csv(sys.argv[1], dtype={'time': str, 'conf': str, 'date': str}, low_memory=False)
        for c in sorted(ex.country.unique()):
            res.append(run(c, prep(ex[ex.country == c])))
    out = {'configuration': 'identical for every country: weight k(calendar month, 0.5-deg cell), '
                            'empirical-Bayes shrinkage alpha=20, fitted on the 2012-2021 MODIS/S-NPP overlap, '
                            'leave-one-year-out cross-validation. No per-country tuning.',
           'note': 'Bangladesh is re-run here under the month-as-season configuration, so its numbers differ '
                   'slightly from the five-season model shipped in the app.',
           'countries': res}
    json.dump(out, open('pipeline/travel_test.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))
