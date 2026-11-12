import sys, json; sys.path.insert(0,'analysis')
import pandas as pd, numpy as np
d=pd.read_pickle('analysis/detections.pkl')
c=json.load(open('analysis/coefficients.json'))
# equivalence weight e per detection (S-NPP-equivalent per detection)
exec(open('analysis/harmonize.py').read().split('# ---- per-detection weights')[0].split('ov=f[f.y.isin')[0])  # defs only
from common import fires
f=d[(d.type==0)&d.confok]
ov=f[f.y.isin(TRAIN)&(f.ym>=SWITCH)]
J_main,ks_main=fit_k(ov,['MA','MT']); J_terra,ks_terra=fit_k(ov,['MT'])
J_n20,ks_n20=fit_k(f[f.y.between(2019,2025)],['V1'])
e=np.ones(len(d),dtype=np.float32)
isM=d.src.isin(['MA','MT']).values; preA=(d.ym<AQUA_START).values
e[isM&~preA]=apply_k(d[isM&~preA],J_main,ks_main)
e[isM&preA]=apply_k(d[isM&preA],J_terra,ks_terra)
v1=(d.src=='V1').values; e[v1]=apply_k(d[v1],J_n20,ks_n20)
d['e']=e
d=d.sort_values(['date','time']).reset_index(drop=True)
n=len(d)
day=((d.dt-pd.Timestamp('2000-01-01')).dt.days).astype(np.uint16).values
tm=d.time.astype(int).astype(np.uint16).values
srcmap={'MT':0,'MA':1,'VS':2,'V1':3}; src=d.src.map(srcmap).astype(np.uint8).values
def cf(row_src,conf):
    if row_src in ('MA','MT'): return int(conf)
    return {'l':201,'n':202,'h':203}[conf]
conf=np.array([cf(s,x) for s,x in zip(d.src,d.conf)],dtype=np.uint8)
typ=d.type.astype(np.uint8).values
flags=((d.daynight=='N').values.astype(np.uint8)|(((d.type==0)&d.confok).values.astype(np.uint8)<<1)).astype(np.uint8)
parts=[d.lat.values.astype(np.float32),d.lon.values.astype(np.float32),d.frp.values.astype(np.float32),d.e.values.astype(np.float32),day,tm,src,conf,typ,flags]
names=['lat','lon','frp','e','day','time','src','conf','type','flags']
off=0; layout={}
buf=bytearray()
for nm,a in zip(names,parts):
    layout[nm]={'offset':len(buf),'dtype':str(a.dtype)}; buf+=a.tobytes()
open('app/points.bin','wb').write(buf)
meta=dict(n=n,layout=layout,epoch='2000-01-01',src_codes=['MODIS Terra','MODIS Aqua','VIIRS S-NPP','VIIRS NOAA-20'],
  switch_ym=201202, aqua_start_ym=200207, first_ym=200011, last_ym=202512,
  seasons=c['seasons'], season_ratio_modis=c['season_ratio_modis'], season_ratio_terra_only=c['season_ratio_terra_only'],
  season_ratio_noaa20=c['season_ratio_noaa20'], uncertainty=c['uncertainty'], validation=c['validation'],
  method=c['method'], filters=c['filters'], data_source=c['data_source'])
import math
def _clean(o):
    if isinstance(o,dict): return {k:_clean(v) for k,v in o.items()}
    if isinstance(o,list): return [_clean(v) for v in o]
    if isinstance(o,float) and (math.isnan(o) or math.isinf(o)): return None
    return o
json.dump(_clean(meta),open('app/meta.json','w'),allow_nan=False)
# reference totals for app verification
d['ym']=d.y*100+d.m
inc=(d.type==0)&d.confok
ser=np.where(inc&(((d.src.isin(['MA','MT']))&(d.ym<SWITCH))|((d.src=='VS')&(d.ym>=SWITCH))),d.e,0)
d['ser']=ser
ref=d.groupby('y').ser.sum().round(1)
json.dump({str(k):float(v) for k,v in ref.items()},open('analysis/ref_annual.json','w'))
print(n,len(buf)/1e6,'MB'); print(ref.to_dict())
