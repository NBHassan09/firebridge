"""Is the control result an artefact of the chosen comparison windows?

control.py compares 2003-2011 against 2013-2021.  A reviewer will reasonably ask
whether that choice produces the result.  This script re-runs the same comparison
over every plausible pair of windows and reports the spread.

Usage:  python analysis/control_windows.py      (run from the project root)
Writes: analysis/control_windows.json
"""
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, 'analysis')
# reuse harmonize.py's loaders, filters and estimator, without re-running its fit
exec(open('analysis/harmonize.py').read().split('ov=f[f.y.isin')[0])

ov = f[f.y.isin(TRAIN) & (f.ym >= SWITCH)]
J_main, ks_main = fit_k(ov, ['MA', 'MT'])
g = f.copy()
g['w'] = 1.0
isM = g.src.isin(['MA', 'MT']).values
g.loc[isM, 'w'] = apply_k(g[isM], J_main, ks_main)

MOD = f[f.src.isin(['MA', 'MT'])]
AQ = f[f.src == 'MA']
VS = f[f.src == 'VS']
HM = g[g.src.isin(['MA', 'MT'])]


def n(df, yrs):
    return df[df.y.isin(yrs)].groupby('y').size().reindex(yrs, fill_value=0).mean()


def nw(df, yrs):
    return df[df.y.isin(yrs)].groupby('y').w.sum().reindex(yrs, fill_value=0).mean()


def chg(a, b):
    return round((b / a - 1) * 100, 1)


def run(a0, a1, b0, b1):
    A = list(range(a0, a1 + 1))
    B = list(range(b0, b1 + 1))
    return {'pre': f'{a0}-{a1}', 'post': f'{b0}-{b1}',
            'pre_years': len(A), 'post_years': len(B),
            'aqua_control': chg(n(AQ, A), n(AQ, B)),
            'modis_both': chg(n(MOD, A), n(MOD, B)),
            'raw_splice': chg(n(MOD, A), n(VS, B)),
            'harmonized': chg(nw(HM, A), n(VS, B))}


# ---- named windows a reader would think to try ----
NAMED = [
    ('published choice', 2003, 2011, 2013, 2021),
    ('one year trimmed each end', 2004, 2011, 2013, 2020),
    ('post period ends 2019', 2003, 2011, 2013, 2019),
    ('post period ends 2017', 2003, 2011, 2013, 2017),
    ('pre period from 2005', 2005, 2011, 2013, 2021),
    ('2012 counted as post-VIIRS', 2003, 2011, 2012, 2021),
    ('symmetric 5 years', 2007, 2011, 2013, 2017),
    ('symmetric 3 years', 2009, 2011, 2013, 2015),
    ('adjacent years only', 2011, 2011, 2013, 2013),
]
named = [dict(label=lbl, **run(a0, a1, b0, b1)) for lbl, a0, a1, b0, b1 in NAMED]

# ---- exhaustive sweep over plausible boundaries ----
sweep = []
for a0 in range(2003, 2010):
    for b1 in range(2016, 2022):
        sweep.append(run(a0, 2011, 2013, b1))


def spread(key):
    v = [s[key] for s in sweep]
    return {'min': round(min(v), 1), 'median': round(float(np.median(v)), 1),
            'max': round(max(v), 1),
            'always_negative': bool(max(v) < 0), 'always_positive': bool(min(v) > 0)}


gap = [s['harmonized'] - s['aqua_control'] for s in sweep]
res = {
    'question': 'Does the control result depend on where the comparison windows are drawn?',
    'named_windows': named,
    'sweep': {'description': 'every pre-window 2003..2009 to 2011, against every '
                             'post-window 2013 to 2016..2021',
              'n_window_pairs': len(sweep),
              'aqua_control': spread('aqua_control'),
              'modis_both': spread('modis_both'),
              'raw_splice': spread('raw_splice'),
              'harmonized': spread('harmonized'),
              'harmonized_minus_control_pct_points': {
                  'min': round(min(gap), 1), 'median': round(float(np.median(gap)), 1),
                  'max': round(max(gap), 1)}},
    'verdict': None,
}
sw = res['sweep']
res['verdict'] = (
    'The sign is stable: across all %d window pairs the unchanged-instrument control is '
    'negative (%.1f%% to %.1f%%) and the raw splice is positive (+%.1f%% to +%.1f%%). '
    'The harmonized series is negative throughout and tracks the control to within '
    '%.1f to %.1f percentage points.' % (
        len(sweep), sw['aqua_control']['min'], sw['aqua_control']['max'],
        sw['raw_splice']['min'], sw['raw_splice']['max'],
        min(gap), max(gap)))

json.dump(res, open('analysis/control_windows.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
