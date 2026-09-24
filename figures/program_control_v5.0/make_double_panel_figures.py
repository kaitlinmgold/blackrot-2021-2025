"""Plot severity (A) and incidence (B) control from the prepared tables.
Year colors and program membership are shared between panels.
"""
from pathlib import Path
import os, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
OUT=Path(__file__).resolve().parent
p=pd.read_csv(OUT/'program_key.csv');m=pd.read_csv(OUT/'product_group_membership.csv')
colors={2021:'#0072B2',2022:'#D55E00',2023:'#009E73',2024:'#AA4499',2025:'#8A6811'}
shapes={'Repeated single product':('o','Single product'),'Repeated tank mix':('s','Tank mix'),'Alternating rotation':('^','Rotation'),'Complex sequence':('D','Complex sequence'),'Formulated hybrid product':('P','Premix')}
contexts={'Single biological product':('o','Single biological product'),'Biological combination / sequence':('^','Biological combination / sequence'),'Includes conventional fungicide':('s','Includes conventional fungicide')}
classes=['Biological only','Mixed','Conventional only']
families=['Howler','Howler EVO','Regalia +\nStargus','LifeGard WG','Oso 5SC','Curezin','OxiDate*','PERpose*']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':18,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
# Verify both endpoints against the existing raw means, without fitting any models.
for ep in ['Severity','Incidence']:
 expected=100*(1-p[f'{ep}_pct']/p[f'UTC_{ep}_pct'])
 assert np.allclose(expected,p[f'{ep}_control_pct'])
 assert np.isfinite(p[f'{ep}_control_pct']).all()
 for _,r in m.iterrows():
  assert np.isclose(r[f'{ep}_control_pct'],p.set_index('Public_Program_ID').loc[r.Public_Program_ID,f'{ep}_control_pct'])
ymin=min(-20,20*np.floor((min(p.Severity_control_pct.min(),p.Incidence_control_pct.min())-3)/20))
positions=[];figures=[]

def draw_group(ax,f,x,endpoint,kind,category):
 col=f'{endpoint}_control_pct'
 if len(f)>=5:
  ax.boxplot([f[col]],positions=[x],widths=.58,whis=(0,100),showfliers=False,patch_artist=True,manage_ticks=False,boxprops={'facecolor':'#EEF1F3','edgecolor':'#89949C','linewidth':1},medianprops={'color':'#46535D','linewidth':1.5},whiskerprops={'color':'#89949C'},capprops={'color':'#89949C'})
 used=[]
 for _,r in f.sort_values([col,'Public_Program_ID']).iterrows():
  v=float(r[col]); step=.055 if kind=='class' else .075
  candidates=[0]+[s*i*step for i in range(1,7 if kind=='class' else 5) for s in [-1,1]]
  valid=[z for z in candidates if all(abs(v-vv)>4 or abs(z-zz)>step*.95 for zz,vv in used)]
  assert valid, 'Insufficient room for non-overlapping dots'
  off=valid[0];used.append((off,v))
  marker=(shapes[r.Architecture_Clean] if kind=='class' else contexts[r.Context])[0]
  ax.scatter(x+off,v,s=100 if kind=='class' else 85,marker=marker,color=colors[int(r.Year)],edgecolors='white',linewidths=.65,zorder=4)
  positions.append({'Figure':kind,'Panel':'A' if endpoint=='Severity' else 'B','Endpoint':endpoint,'Category':category.replace('\n',' '),'Program_ID':r.Public_Program_ID,'Year':int(r.Year),'Control_pct':v,'X_position':x+off,'Marker':marker})

def make(kind):
 isclass=kind=='class'; frame=p if isclass else m; cats=classes if isclass else families; key='Program_Class_Clean' if isclass else 'Product_group'
 fig,axs=plt.subplots(2,1,figsize=(14,12),sharex=True,sharey=True)
 fig.subplots_adjust(left=.115,right=.975,bottom=.115,top=.77,hspace=.25)
 fig.suptitle('Control across spray strategies' if isclass else 'Control of programs containing recurring products',x=.055,y=.976,ha='left',fontsize=24,fontweight='bold',color='#20313B')
 
 handles=[Line2D([],[],marker='o',linestyle='none',markerfacecolor=c,markeredgecolor='white',markersize=11,label=str(y)) for y,c in colors.items()]
 fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.53,.94),ncol=5,frameon=False,fontsize=18)
 shapehandles=[Line2D([],[],marker=mk,linestyle='none',markerfacecolor='#58656E',markeredgecolor='white',markersize=11,label=lab) for mk,lab in (shapes if isclass else contexts).values()]
 fig.legend(handles=shapehandles,loc='upper center',bbox_to_anchor=(.53,.892),ncol=3 if isclass else 2,frameon=False,fontsize=17, columnspacing=1.4, handletextpad=.35)
 for ax,letter,ep in zip(axs,['A','B'],['Severity','Incidence']):
  for x,c in enumerate(cats):draw_group(ax,frame[frame[key]==c],x,ep,kind,c)
  ax.set_title(f'{letter}   {ep} control',loc='left',fontsize=22,fontweight='bold',pad=10,color='#20313B')
  ax.set_ylim(ymin,105);ax.set_yticks(np.arange(ymin,101,20));ax.set_ylabel(f'{ep} control (%)',fontsize=18)
  ax.set_axisbelow(True);ax.grid(axis='y',color='#E4E8EB',linewidth=.7);ax.axhline(0,color='#59656E',lw=.9,ls='--')
  for side in ['left','bottom']:ax.spines[side].set_color('#8E989F')
  ax.set_xlim(-.55,len(cats)-.45)
 labels=['Biological-only','Mixed','Conventional-only'] if isclass else cats
 axs[-1].set_xticks(range(len(cats)),[f'{lab}\n(n = {sum(frame[key]==c)})' for lab,c in zip(labels,cats)],fontsize=18 if isclass else 16)
 notes=['Boxes show the middle 50%, median, and full observed range'+('.' if isclass else '; boxes appear for n ≥ 5 program-years.'),
        'Control (%) = 100 × (1 − program mean / same-year untreated mean), calculated separately for severity and incidence.']
 if isclass:
  notes+=['Biological-only includes metal products under the manuscript’s classification; mixed programs include biological and conventional inputs.',
          '59 identifiable treated program-years per panel; two uncertain 2022 identities omitted. Descriptive summaries; no new statistical tests.']
 else:
  notes+=['A program may appear in several product groups. Dots describe whole-program performance; details are in the program key.',
          '*OxiDate and PERpose group recorded name variants. Howler and Howler EVO are separate; Ensendo is not included under Howler.']
 # Explanatory notes are preserved in the manuscript captions and accompanying caption file.
 stem='Figure_4_control_by_class_v5.0' if isclass else 'Figure_5_control_by_product_v5.0'
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{stem}.{ext}',dpi=220,facecolor='white')
 figures.append(fig)
 return stem
stems=[make('class'),make('product')]
q=pd.DataFrame(positions)
assert q.groupby(['Figure','Panel']).size().to_dict()=={('class','A'):59,('class','B'):59,('product','A'):42,('product','B'):42}
for kind in ['class','product']:
 a=q[(q.Figure==kind)&(q.Panel=='A')];b=q[(q.Figure==kind)&(q.Panel=='B')]
 assert set(zip(a.Program_ID,a.Category,a.Marker))==set(zip(b.Program_ID,b.Category,b.Marker))
q.to_csv(OUT/'double_panel_plotted_values.csv',index=False)
with PdfPages(OUT/'BlackRot_double_panel_figures.pdf') as pdf:
 for fig in figures:pdf.savefig(fig)
(OUT/'double_panel_validation.json').write_text(json.dumps({'class_program_years_per_panel':59,'product_memberships_per_panel':42,'product_unique_program_years':m.Public_Program_ID.nunique(),'panels':['A: severity control','B: incidence control'],'both_endpoints_verified_against_existing_raw_means':True,'identical_program_membership_and_marker_mapping_between_panels':True,'common_y_limits':[float(ymin),105],'no_new_statistical_tests':True},indent=2))
print('Both double-panel figures saved; values and membership verified.')
print(p[['Severity_control_pct','Incidence_control_pct']].agg(['min','max']).to_string())
