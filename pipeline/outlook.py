"""Fire-season outlook: predict Mar-Apr hotspots from conditions known by 1 March.
Leave-one-year-out over 2003-2025 (MODIS-estimated + VIIRS years)."""
import json, numpy as np, pandas as pd
from scipy.stats import spearmanr
d=pd.read_pickle('analysis/detections.pkl'); c=json.load(open('app/context.json'))
inc=(d.type==0)&d.confok
d['ser']=np.where(inc&(((d.src.isin(['MA','MT']))&(d.ym<201202))|((d.src=='VS')&(d.ym>=201202))),d.w,0)
def clim(k,box):
    pts=c['params'][k]['points']
    if box: pts=[p for p in pts if box[0]-.3<=p['lat']<=box[1]+.3 and box[2]-.35<=p['lon']<=box[3]+.35]
    return np.nanmean(np.array([[np.nan if x is None else x for x in p['v']] for p in pts],dtype=float),axis=0)
ix=lambda y,m:(y-2001)*12+m-1
DIM=lambda y,m:pd.Period(f'{y}-{m:02d}').days_in_month
def table(box):
    dd=d if box is None else d[d.lat.between(box[0],box[1])&d.lon.between(box[2],box[3])]
    Hm=dd.groupby('ym').ser.sum(); pr=clim('PRECTOTCORR',box); gw=clim('GWETTOP',box); rh=clim('RH2M',box); tx=clim('T2M_MAX',box)
    rows=[]
    for y in range(2003,2026):
        rows.append(dict(y=y,fire=Hm.get(y*100+3,0)+Hm.get(y*100+4,0),
          rain_NDJF=sum(pr[ix(yy,m)]*DIM(yy,m) for yy,m in [(y-1,11),(y-1,12),(y,1),(y,2)]),
          rain_JF=sum(pr[ix(y,m)]*DIM(y,m) for m in [1,2]),
          soil_Feb=gw[ix(y,2)], rh_Feb=rh[ix(y,2)], tmax_Feb=tx[ix(y,2)]))
    return pd.DataFrame(rows)
def loyo(R,cols):
    pred=[]; 
    for k in range(len(R)):
        tr=R.drop(R.index[k]); X=np.column_stack([np.ones(len(tr))]+[(tr[c]-tr[c].mean())/tr[c].std() for c in cols]); yv=np.log(tr.fire+1)
        b,*_=np.linalg.lstsq(X,yv,rcond=None); te=R.iloc[k]; x=np.array([1]+[(te[c]-tr[c].mean())/tr[c].std() for c in cols]); pred.append(np.exp(x@b)-1)
    R=R.copy(); R['pred']=pred
    # tercile skill: category of actual vs predicted relative to training terciles (approx: full-sample terciles)
    q=R.fire.quantile([1/3,2/3]).values; cat=lambda v:0 if v<q[0] else (2 if v>q[1] else 1)
    hit=np.mean([cat(a)==cat(p) for a,p in zip(R.fire,R.pred)])
    rho=spearmanr(R.pred,R.fire).statistic; mape=np.median(abs(R.pred-R.fire)/R.fire)
    return rho,hit,mape,R
for name,box in [('Bangladesh',None),('CHT',(21.15,23.75,91.85,92.7)),('Plains west (Barind)',(24.3,25.3,88.1,89.2))]:
    R=table(box); print(name)
    for col in ['rain_NDJF','rain_JF','soil_Feb','rh_Feb','tmax_Feb']:
        print(f'   {col:10s} rho={spearmanr(R[col],R.fire).statistic:+.2f}')
    for cols in [['soil_Feb'],['rain_NDJF'],['rh_Feb'],['soil_Feb','rh_Feb'],['rain_NDJF','rh_Feb'],['soil_Feb','tmax_Feb']]:
        rho,hit,mape,_=loyo(R,cols); print(f'   LOYO {"+".join(cols):20s} rho={rho:.2f} tercile hit={hit:.2f} (chance 0.33) median err={mape:.2f}')
