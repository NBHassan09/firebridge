"""MODIS -> VIIRS(S-NPP) harmonization of NASA FIRMS hotspots, Bangladesh.
Reference series: VIIRS S-NPP 375 m (VNP14IMG, standard processing), Feb 2012 onward.
MODIS (MCD14ML C6.1) detections are converted to S-NPP-equivalent counts with a
per-detection weight k(season, 0.5-degree cell), fitted on 2012-2021 overlap with
empirical-Bayes shrinkage toward the season-wide ratio.  NOAA-20 VIIRS is linked to
S-NPP the same way so the series continues after S-NPP ended (2 Nov 2026)."""
import sys, json; sys.path.insert(0,'analysis')
from common import *
from scipy.stats import spearmanr
ALPHA=20; RES=0.5
TRAIN=list(range(2012,2022)); SWITCH=201202   # first month of S-NPP-based series
AQUA_START=200207
SEAS={1:'NovFeb',2:'NovFeb',11:'NovFeb',12:'NovFeb',3:'Mar',4:'Apr',5:'May'}
d=load(); d['s']=d.m.map(SEAS).fillna('JunOct')
d['cell']=(np.floor(d.lat/RES)).astype(int)*1000+(np.floor(d.lon/RES)).astype(int)
F=(d.type==0)&d.confok
f=d[F]

def fit_k(train, pred, ref='VS'):
    X=train[train.src.isin(pred)].groupby(['s','cell']).size().rename('X')
    Y=train[train.src==ref].groupby(['s','cell']).size().rename('Y')
    J=pd.concat([X,Y],axis=1).fillna(0)
    ks=(J.groupby(level=0).Y.sum()/J.groupby(level=0).X.sum())
    J['k0']=ks.reindex(J.index.get_level_values(0)).values
    J['k']=(J.Y+ALPHA*J.k0)/(J.X+ALPHA)
    fac=J.groupby(level=0).Y.sum()/(J.k*J.X).groupby(level=0).sum()
    J['k']=J.k*fac.reindex(J.index.get_level_values(0)).values
    return J, ks
def apply_k(df, J, ks):
    kd=J.k.to_dict()
    return np.array([kd.get((s,c),ks[s]) for s,c in zip(df.s,df.cell)])

ov=f[f.y.isin(TRAIN)&(f.ym>=SWITCH)]
# ---- CV with residual collection (leave-one-year-out) ----
def cv(pred, years=TRAIN, data=ov):
    rows=[]
    for yt in years:
        J,ks=fit_k(data[data.y!=yt],pred); te=data[data.y==yt]
        xm=te[te.src.isin(pred)].copy(); xm['w']=apply_k(xm,J,ks)
        P=xm.groupby(['m','cell']).w.sum(); O=te[te.src=='VS'].groupby(['m','cell']).size()
        R=pd.concat([P.rename('P'),O.rename('O')],axis=1).fillna(0); R['y']=yt; rows.append(R)
    return pd.concat(rows)
def metrics(R):
    nat=R.groupby(['y','m']).sum(); nb=nat[nat.O>=50]
    ann=R.groupby(['y','cell']).sum(); ann=ann[ann.O>=30]
    return dict(national_month_median_abs_pct_err=round(float((abs(nb.P-nb.O)/nb.O).median()),3),
                cell_month_WAPE=round(float(abs(R.P-R.O).sum()/R.O.sum()),3),
                cell_annual_spearman=round(float(spearmanr(ann.P,ann.O).statistic),3),
                total_bias=round(float(R.P.sum()/R.O.sum()-1),3))
R_main=cv(['MA','MT']); R_terra=cv(['MT'])
# baseline comparisons for the write-up
def cv_simple(pred,mode):
    rows=[]
    for yt in TRAIN:
        tr=ov[ov.y!=yt]; te=ov[ov.y==yt]
        X=tr[tr.src.isin(pred)]; Y=tr[tr.src=='VS']
        if mode=='global': kk=lambda df: np.full(len(df),len(Y)/len(X))
        else:
            ks=Y.groupby('s').size()/X.groupby('s').size(); kk=lambda df: df.s.map(ks).values
        xm=te[te.src.isin(pred)].copy(); xm['w']=kk(xm)
        P=xm.groupby(['m','cell']).w.sum(); O=te[te.src=='VS'].groupby(['m','cell']).size()
        Rr=pd.concat([P.rename('P'),O.rename('O')],axis=1).fillna(0); Rr['y']=yt; rows.append(Rr)
    return metrics(pd.concat(rows))
compare={'single_global_ratio':cv_simple(['MA','MT'],'global'),
         'seasonal_ratio':cv_simple(['MA','MT'],'season'),
         'seasonal_x_cell_shrinkage (chosen)':metrics(R_main),
         'terra_only_seasonal_x_cell':metrics(R_terra)}
print(json.dumps(compare,indent=1))

# ---- uncertainty: empirical 10th/90th percentiles of observed/harmonized, binned by harmonized count,
# pooled over CV residuals at several aggregation levels (cell, 2x2, 4x4 cells, national; month and year)
BINS=[0,3,10,30,100,300,1000,1e9]
def levels(R):
    Rg=R.reset_index(); out=[]
    for g in [1,2,4,100]:
        Rg['g']=(Rg.cell//1000//g)*1000+(Rg.cell%1000)//g
        out.append(Rg.groupby(['y','m','g'])[['P','O']].sum())
        out.append(Rg.groupby(['y','g'])[['P','O']].sum())
    A=pd.concat(out); return A[A.P>0.5]
def unc(R):
    A=levels(R); A['b']=pd.cut(A.P,BINS,right=False,labels=False); q=[]
    for b in range(len(BINS)-1):
        x=A[A.b==b]; 
        if len(x)<15: q.append(q[-1] if q else [0,3]); continue
        r=x.O/x.P; q.append([round(float(r.quantile(.1)),3),round(float(r.quantile(.9)),3)])
    lo=A.P*np.array([q[b][0] for b in A.b]); hi=A.P*np.array([q[b][1] for b in A.b])
    cover=float(((A.O>=lo)&(A.O<=hi)).mean())
    return dict(bins=BINS[:-1],q10_q90_ratio=q,coverage_80=round(cover,3),n=int(len(A)))
U_main=unc(R_main); U_terra=unc(R_terra); print('unc',U_main,U_terra)

# ---- final fits on all training years
J_main,ks_main=fit_k(ov,['MA','MT']); J_terra,ks_terra=fit_k(ov,['MT'])
# NOAA-20 -> S-NPP link (full overlap years 2019-2025)
ovn=f[f.y.between(2019,2025)]
J_n20,ks_n20=fit_k(ovn,['V1'])
# NOAA-20 linkage validation: LOYO
Rn=[]
for yt in range(2019,2026):
    J,ks=fit_k(ovn[ovn.y!=yt],['V1']); te=ovn[ovn.y==yt]; xm=te[te.src=='V1'].copy(); xm['w']=apply_k(xm,J,ks)
    P=xm.groupby(['m','cell']).w.sum(); O=te[te.src=='VS'].groupby(['m','cell']).size()
    r=pd.concat([P.rename('P'),O.rename('O')],axis=1).fillna(0); r['y']=yt; Rn.append(r)
Rn=pd.concat(Rn); n20_metrics=metrics(Rn); print('n20',n20_metrics)

# ---- out-of-sample drift test: model trained 2012-2021 applied to 2022-2025 MODIS
post=f[f.y.between(2022,2025)].copy()
xm=post[post.src.isin(['MA','MT'])].copy(); xm['w']=apply_k(xm,J_main,ks_main)
drift=pd.DataFrame({'modis_based':xm.groupby('y').w.sum().round(0),'snpp_observed':post[post.src=='VS'].groupby('y').size()})
drift['ratio']=(drift.modis_based/drift.snpp_observed).round(2)
med_time=d[(d.daynight=='D')&F&d.src.isin(['MA','MT'])].groupby(['y','src']).time.apply(lambda t:(t.str[:2].astype(int)+t.str[2:].astype(int)/60+6).median()).unstack().round(2)
print(drift); print(med_time.loc[2019:])

# ---- per-detection weights for the harmonized series
w=np.zeros(len(d))
m=F.values
isMODIS=d.src.isin(['MA','MT']).values
pre_aqua=(d.ym<AQUA_START).values; pre_switch=(d.ym<SWITCH).values
sel=m&isMODIS&pre_switch&~pre_aqua
w[sel]=apply_k(d[sel],J_main,ks_main)
sel=m&(d.src=='MT').values&pre_aqua
w[sel]=apply_k(d[sel],J_terra,ks_terra)
sel=m&(d.src=='VS').values&~pre_switch
w[sel]=1.0
d['w']=w
d['era']=np.select([d.ym<AQUA_START,d.ym<SWITCH],[0,1],2)   # 0 Terra-only, 1 MODIS T+A, 2 VIIRS S-NPP
# sanity: monthly harmonized national series
H=d.groupby('ym').w.sum()
print('annual harmonized', H.groupby(H.index//100).sum().round(0).to_dict())

coef=dict(
  method="Per-detection weights convert MODIS (Terra+Aqua, C6.1) active-fire detections into VIIRS S-NPP 375 m-equivalent counts. Weight k = ratio of S-NPP to MODIS detections in the same season and 0.5-degree cell during the 2012-2021 overlap, shrunk toward the season-wide ratio (alpha=20 MODIS detections) and rescaled to preserve season totals. Before Aqua (to Jun 2002) a Terra-only ratio is used. From Feb 2012 the series is S-NPP itself. NOAA-20 is linked to S-NPP with the same estimator (2019-2025) for continuity after S-NPP ended on 2 Nov 2026.",
  filters="type=0 (presumed vegetation fire) only; MODIS confidence>=30; VIIRS confidence nominal or high. Static land sources (type 2) and offshore (type 3) are excluded from counts but kept in exports.",
  seasons={"NovFeb":[11,12,1,2],"Mar":[3],"Apr":[4],"May":[5],"JunOct":[6,7,8,9,10]},
  season_ratio_modis={k:round(float(v),3) for k,v in ks_main.items()},
  season_ratio_terra_only={k:round(float(v),3) for k,v in ks_terra.items()},
  season_ratio_noaa20={k:round(float(v),3) for k,v in ks_n20.items()},
  cell_ratios_modis={f"{s}|{c}":round(float(v),3) for (s,c),v in J_main.k.items()},
  uncertainty={"modis_era":U_main,"terra_only_era":U_terra,"form":"80% range = [H*q10, H*q90] using the bin of H"},
  validation=dict(cv="leave-one-year-out over 2012-2021",model_comparison=compare,noaa20_link=n20_metrics,
                  drift_test={str(k):v for k,v in drift.to_dict(orient='index').items()},
                  modis_median_overpass_local={str(y):r for y,r in med_time.to_dict(orient='index').items()}),
  eras={"0":"Terra-only (Nov 2000-Jun 2002)","1":"MODIS Terra+Aqua harmonized (Jul 2002-Jan 2012)","2":"VIIRS S-NPP observed (Feb 2012-Dec 2025)"},
  data_source="NASA FIRMS country yearly summaries (standard processing): MODIS C6.1 MCD14ML, VIIRS S-NPP VNP14IMG, VIIRS NOAA-20 VJ114IMG; firms.modaps.eosdis.nasa.gov/country/")
json.dump(coef,open('analysis/coefficients.json','w'),indent=1)
d.to_pickle('analysis/detections.pkl')
