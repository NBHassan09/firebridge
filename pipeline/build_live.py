import sys, json, math, numpy as np, pandas as pd
sys.path.insert(0,'pipeline'); sys.path.insert(0,'analysis')
from live_common import Processor, season
# NOAA-21 link: fit on 2026 overlap with S-NPP (Jan-Oct, Bangladesh only) with the same estimator
n=pd.read_csv('data_in/nrt.csv.gz',dtype={'date':str,'time':str,'conf':str}).drop_duplicates(['src','date','time','lat','lon'])
n=n[~n.date.str.startswith('2025')]
json.dump({'n20':json.load(open('live/link_n20.json'))['n20'],'n21':{'season':{},'cell':{}}},open('live/link_coeffs.json','w'))
P=Processor('live/live_ref.npz','live/link_coeffs.json','live/landcover_2024.npy.gz')
bd,_=P.in_bd(n.lat.values,n.lon.values); n=n[bd].copy()
n['m']=n.date.str[4:6].astype(int); n['ym']=n.date.str[:6].astype(int)
n['type']=[2 if P.static_cell(a,b) in P.static else 0 for a,b in zip(n.lat,n.lon)]
f=n[(n.type==0)&(n.conf!='l')&(n.ym<=202610)].copy(); f['s']=f.m.map(season); f['cell']=(np.floor(f.lat/0.5)).astype(int)*1000+(np.floor(f.lon/0.5)).astype(int)
def fit(pred):
    X=f[f.src==pred].groupby(['s','cell']).size().rename('X'); Y=f[f.src=='VS'].groupby(['s','cell']).size().rename('Y')
    J=pd.concat([X,Y],axis=1).fillna(0); ks=J.groupby(level=0).Y.sum()/J.groupby(level=0).X.sum()
    J['k0']=ks.reindex(J.index.get_level_values(0)).values; J['k']=(J.Y+20*J.k0)/(J.X+20)
    fac=J.groupby(level=0).Y.sum()/(J.k*J.X).groupby(level=0).sum(); J['k']=J.k*fac.reindex(J.index.get_level_values(0)).values
    return {'season':{k:round(float(v),4) for k,v in ks.items()},'cell':{f'{s}|{c}':round(float(v),4) for (s,c),v in J.k.items()}}
n21=fit('V2'); n20_nrt=fit('V1')
co=json.load(open('live/link_coeffs.json')); co['n21']=n21; json.dump(co,open('live/link_coeffs.json','w'))
print('NOAA-21 season ratios (S-NPP per NOAA-21):',n21['season']); print('NOAA-20 in 2026 NRT:',n20_nrt['season'],' vs 2019-2025 SP link:',co['n20']['season'])
# monthly check: S-NPP vs NOAA-20 converted, 2026
P=Processor('live/live_ref.npz','live/link_coeffs.json','live/landcover_2024.npy.gz')
rows=[[r.src,r.date,r.time,r.lat,r.lon,r.conf,r.frp,r.daynight] for r in n.itertuples()]
out=P.process(rows)
L=pd.DataFrame(out,columns=['src','date','time','lat','lon','conf','frp','dn','type','e','lc']); L['ym']=L.date.str[:6]
inc=L[(L.type==0)&(L.conf!='l')]
tab=inc.groupby(['ym','src']).e.sum().unstack().round(0); tab['n20/snpp']=(tab.V1/tab.VS).round(2); tab['n21/snpp']=(tab.V2/tab.VS).round(2)
print(tab.to_string())
print('static flagged',(L.type==2).sum(),'of',len(L))
L=L.drop(columns='ym')
json.dump({'fields':['src','date','time','lat','lon','conf','frp','daynight','type','snpp_equiv_weight','igbp'],
   'source':'NASA GIBS VIIRS_{SNPP,NOAA20,NOAA21}_Thermal_Anomalies_375m_All vector tiles (FIRMS NRT), assigned to Bangladesh with FIRMS country reference points',
   'rows':L.values.tolist(),'updated':'2026-11-12','last_day':L.date.max()},open('live/live.json','w'),separators=(',',':'))
# POWER 2026
pw=json.load(open('data_in/power2026.json')); json.dump({'source':pw['source'],'params':pw['params']},open('live/live_power.json','w'),separators=(',',':'))
print('rows',len(L),'last',L.date.max())
