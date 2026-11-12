import pandas as pd, numpy as np
def load():
    d=pd.read_csv('data_in/firms_bd.csv.gz',dtype={'time':str,'conf':str,'date':str})
    d['dt']=pd.to_datetime(d.date,format='%Y%m%d'); d['y']=d.dt.dt.year; d['m']=d.dt.dt.month
    isM=d.src.isin(['MT','MA'])
    d['confok']=np.where(isM, pd.to_numeric(d.conf,errors='coerce')>=30, d.conf!='l')
    d['ym']=d.y*100+d.m
    return d
def fires(d): return d[(d.type==0)&d.confok].copy()
