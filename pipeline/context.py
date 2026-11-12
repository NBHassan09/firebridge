"""Attach NASA context to the FIRMS detections.
1) MODIS land cover (MCD12Q1 IGBP, via GIBS) for each detection, using the land cover of the detection's year (clamped to 2001-2024).
2) NASA POWER monthly climate per grid point, trimmed to points over Bangladesh, written to app/context.json."""
import json, numpy as np, pandas as pd
ctx=json.load(open('data_in/nasa_context.json')); LC=np.load('data_in/landcover.npy')
lcm=ctx['landcover']; Y0=lcm['years'][0]; NY=len(lcm['years']); H,W=lcm['height'],lcm['width']; bb=lcm['bbox']
GROUPS={'crop':[12,14],'forest':[1,2,3,4,5],'shrub':[6,7,8,9,10],'built':[13],'wet':[11,17],'other':[0,15,16,254,255]}
def igbp(lat,lon,year):
    yi=np.clip(year-Y0,0,NY-1)
    r=np.clip(((bb['north']-lat)/(bb['north']-bb['south'])*H).astype(int),0,H-1)
    c=np.clip(((lon-bb['west'])/(bb['east']-bb['west'])*W).astype(int),0,W-1)
    return LC[yi,r,c]
def power_json(det):
    P=ctx['power']; out={'source':P['source'],'first':'200101','last':'202512','params':{}}
    lat=det.lat.values; lon=det.lon.values
    for k,v in P['params'].items():
        pts=[]
        half_lat,half_lon=(0.5,0.5) if k=='AOD_55' else (0.25,0.3125)
        for q in v['points']:
            n=int(((abs(lat-q['lat'])<=half_lat)&(abs(lon-q['lon'])<=half_lon)).sum())
            if n==0: continue
            pts.append({'lat':q['lat'],'lon':q['lon'],'n':n,'v':[None if x is None else round(x,3) for x in q['v']]})
        out['params'][k]={'units':v['units'],'name':v['longname'],'half':[half_lat,half_lon],'points':pts}
    return out
if __name__=='__main__':
    d=pd.read_pickle('analysis/detections.pkl')
    d['igbp']=igbp(d.lat.values,d.lon.values,d.y.values)
    inv={c:g for g,cs in GROUPS.items() for c in cs}; d['lcg']=d.igbp.map(inv)
    f=d[(d.type==0)&d.confok]
    print('share of VIIRS S-NPP fire detections by land cover group, 2012-2025'); print((f[f.src=='VS'].lcg.value_counts(normalize=True)*100).round(1))
    print('IGBP classes of S-NPP fires'); print(f[f.src=='VS'].igbp.value_counts().head(10))
    pj=power_json(d); print({k:len(v['points']) for k,v in pj['params'].items()})
    json.dump(pj,open('app/context.json','w'),separators=(',',':'),allow_nan=False)
    d[['igbp']].to_pickle('analysis/igbp.pkl')
