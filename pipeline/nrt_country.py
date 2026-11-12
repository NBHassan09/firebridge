"""Assign GIBS NRT detections to Bangladesh using FIRMS's own country assignment as reference.
Labels: BD = FIRMS Bangladesh country-file detections (VIIRS 2012-2025); foreign = March 2025 GIBS detections
in the bounding box that are not in the Bangladesh country file. Tested on held-out March 2025 days."""
import numpy as np, pandas as pd, math
from scipy.spatial import cKDTree
COS=math.cos(math.radians(23.7))
d=pd.read_pickle('analysis/detections.pkl')
n=pd.read_csv('data_in/nrt.csv.gz',dtype={'date':str,'time':str,'conf':str}).drop_duplicates(['src','date','time','lat','lon'])
m25=n[n.date.str.startswith('202503')].copy()
bd=d[d.src.isin(['VS','V1'])]
# truth for March 2025: is a GIBS point in the BD country file (same sensor, same day, within ~200 m)?
truth=np.zeros(len(m25),bool)
for (src,day),g in m25.groupby(['src','date']):
    ref=bd[(bd.src==src)&(bd.date==day)]
    if len(ref)==0: continue
    t=cKDTree(np.column_stack([ref.lat,ref.lon*COS])); dist,_=t.query(np.column_stack([g.lat,g.lon*COS]))
    truth[m25.index.get_indexer(g.index)]=dist<0.002
m25['bd']=truth
sp=bd[bd.date.str.startswith('202503')]
print('March 2025: GIBS points in box',len(m25),'matched to BD file',truth.sum(),'| BD country file points (VS+V1)',len(sp))
print('  matched share of BD file: %.3f'%(truth.sum()/len(sp)))
days=sorted(m25.date.unique()); train_days=set(days[::2]); test=m25[~m25.date.isin(train_days)]; trainf=m25[m25.date.isin(train_days)&~m25.bd]
hist=bd[~bd.date.str.startswith('202503')]
def build(foreign):
    lat=np.r_[hist.lat.values,foreign.lat.values]; lon=np.r_[hist.lon.values,foreign.lon.values]; lab=np.r_[np.ones(len(hist),bool),np.zeros(len(foreign),bool)]
    return cKDTree(np.column_stack([lat,lon*COS])),lab
for k in [1,3,5,9]:
    tree,lab=build(trainf); _,i=tree.query(np.column_stack([test.lat,test.lon*COS]),k=k)
    pred=(lab[i].reshape(len(test),-1).mean(axis=1)>0.5)
    tp=(pred&test.bd).sum(); fp=(pred&~test.bd).sum(); fn=(~pred&test.bd).sum()
    print(f'k={k}: precision {tp/(tp+fp):.4f} recall {tp/(tp+fn):.4f} | BD count predicted {pred.sum()} vs true {test.bd.sum()} ({(pred.sum()/test.bd.sum()-1)*100:+.2f}%)')
# final reference: all BD history + all March-2025 foreign points
foreign=m25[~m25.bd]
st=d[d.type==2]; cells=set()
for la,lo in zip(st.lat,st.lon):
    r,c=int(round(la/0.004)),int(round(lo/0.004))
    for dr in (-1,0,1):
        for dc in (-1,0,1): cells.add((r+dr,c+dc))
np.savez_compressed('live/live_ref.npz',lat=np.r_[bd.lat.values,foreign.lat.values].astype(np.float32),lon=np.r_[bd.lon.values,foreign.lon.values].astype(np.float32),
    bd=np.r_[np.ones(len(bd),bool),np.zeros(len(foreign),bool)],static_cells=np.array(sorted(cells),dtype=np.int32))
print('reference points: BD',len(bd),'foreign',len(foreign))
