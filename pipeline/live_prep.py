"""Prepare reference files for the NRT pipeline: static-source cells, 2024 land cover, NOAA-20 link coefficients."""
import sys, json, gzip, io, numpy as np, pandas as pd
sys.path.insert(0,'analysis')
exec(open('analysis/harmonize.py').read().split('ov=f[f.y.isin')[0])  # load d, fit_k, constants
f=d[(d.type==0)&d.confok]
J_n20,ks_n20=fit_k(f[f.y.between(2019,2025)],['V1'])
link={'n20':{'season':{k:round(float(v),4) for k,v in ks_n20.items()},'cell':{f'{s}|{c}':round(float(v),4) for (s,c),v in J_n20.k.items()}}}
json.dump(link,open('live/link_n20.json','w'))
# static cells: historic type-2 detections, rounded to 0.004 deg, dilated by one cell
st=d[d.type==2]; cells=set()
for la,lo in zip(st.lat,st.lon):
    r,c=int(round(la/0.004)),int(round(lo/0.004))
    for dr in (-1,0,1):
        for dc in (-1,0,1): cells.add((r+dr,c+dc))
np.save('live/static_cells.npy',np.array(sorted(cells),dtype=np.int32))
lc=np.load('data_in/landcover.npy')[-1]  # 2024
buf=io.BytesIO(); np.save(buf,lc); open('live/landcover_2024.npy.gz','wb').write(gzip.compress(buf.getvalue(),9))
print('n20 seasons',link['n20']['season'],'static cells',len(cells))
