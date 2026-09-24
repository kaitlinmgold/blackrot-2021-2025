"""Derive Figure 4–5 input tables from plot means and application sequences.
Each program-year receives equal weight; no statistical models are fitted.
"""
from pathlib import Path
import hashlib
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT
OUT=Path(__file__).resolve().parent
# Read only: original endpoints, classifications, and application records.
files=[SRC/'data/clean_plot_data_public.csv',SRC/'data/retained_program_codebook_public.csv']
hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
d=pd.read_csv(files[0]);cb=pd.read_csv(files[1]).fillna('')
a=d.groupby(['Year','Public_Program_ID'],as_index=False).agg(N_plots=('Public_Plot_ID','size'),Severity_pct=('Cluster_Severity_Pct','mean'),Incidence_pct=('Cluster_Incidence_Pct','mean'),Clusters_rated=('N_Clusters_Rated','sum'))
a=a.merge(cb,on=['Year','Public_Program_ID'],validate='one_to_one')
assert len(a)==66 and len(d)==264
utc=a[a.Program_Class_Clean=='Untreated'].set_index('Year')
assert len(utc)==5
for ep in ['Severity','Incidence']:
 a[f'UTC_{ep}_pct']=a.Year.map(utc[f'{ep}_pct'])
 a[f'{ep}_control_pct']=np.where(a[f'UTC_{ep}_pct']!=0,100*(1-a[f'{ep}_pct']/a[f'UTC_{ep}_pct']),np.nan)
a['Identity_uncertain']=a.Public_Program_ID.isin(['Y2022_P02','Y2022_P03'])
a['Display_ID']=a.Public_Program_ID.str.replace('Y','',regex=False)
a['Program_number']=a.Public_Program_ID.str.extract(r'_P(\d+)$')[0]
a['Application_sequence']=a.apply(lambda r:'; '.join(f'{i}: {r[f"Application_{i}"]}' for i in range(1,6) if r[f'Application_{i}']),axis=1)
a['Rate_confirmation_pending']=a.Public_Program_ID.isin(['Y2022_P16','Y2023_P07'])
a['Displayed']=~a.Identity_uncertain & (a.Program_Class_Clean!='Untreated')
p=a[a.Displayed].copy();assert len(p)==59

def strategy(r):
 c=r.Program_Class_Clean; ar=r.Architecture_Clean
 if c=='Conventional only':return 'Conventional\nprograms'
 pre='Biological-only' if c=='Biological only' else 'Mixed'
 sub={'Repeated single product':'single product','Repeated tank mix':'tank mix','Alternating rotation':'rotation','Complex sequence':'complex sequence','Formulated hybrid product':'premix'}[ar]
 return pre+'\n'+sub
p['Strategy']=p.apply(strategy,axis=1)
order=['Biological-only\nsingle product','Biological-only\ntank mix','Biological-only\nrotation','Biological-only\ncomplex sequence','Mixed\nrotation','Mixed\ntank mix','Mixed\npremix','Mixed\ncomplex sequence','Conventional\nprograms']
p['Context']=np.where(p.Program_Class_Clean=='Mixed','Includes conventional fungicide',np.where(p.Architecture_Clean=='Repeated single product','Single biological product','Biological combination / sequence'))
# Family membership is based on explicitly recorded applications, not efficacy.
families={'Howler':lambda t:'howler' in t and 'howler evo' not in t,
 'Howler EVO':lambda t:'howler evo' in t,
 'Regalia +\nStargus':lambda t:'regalia' in t and 'stargus' in t,
 'LifeGard WG':lambda t:'lifegard' in t,
 'Oso 5SC':lambda t:'oso 5sc' in t,
 'Curezin':lambda t:'curezin' in t,
 'OxiDate*':lambda t:'oxidate' in t,
 'PERpose*':lambda t:'perpose' in t}
mem=[]
for _,r in p.iterrows():
 for fam,fn in families.items():
  if fn(r.Application_sequence.lower()):
   rec=r.to_dict();rec['Product_group']=fam;mem.append(rec)
m=pd.DataFrame(mem)
# Howler-containing premix Ensendo remains separate; it is not recorded as a Howler application.
a=a.merge(p[['Public_Program_ID','Strategy']],on='Public_Program_ID',how='left')
a.to_csv(OUT/'all_program_control_values.csv',index=False)
p.to_csv(OUT/'program_key.csv',index=False)
m.to_csv(OUT/'product_group_membership.csv',index=False)
utc.reset_index()[['Year','Public_Program_ID','N_plots','Severity_pct','Incidence_pct','Clusters_rated']].to_csv(OUT/'annual_untreated_reference.csv',index=False)
