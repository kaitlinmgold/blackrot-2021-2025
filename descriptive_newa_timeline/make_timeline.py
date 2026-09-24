from pathlib import Path
import json,csv,shutil
from datetime import datetime, timedelta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
out=Path(__file__).resolve().parent
methods=list(csv.DictReader(open(out/'annual_methods_weather.csv')))
rows=[]
for y in range(2021,2026):
 t=json.loads((out/'source_logs'/f'NEWA_Crittenden_{y}_tables.json').read_text())['tables'][3]
 for r in t[1:]:
  if r[6]!='Infection':continue
  s,e=[datetime.strptime(v,'%m/%d/%Y %I:%M %p') for v in r[:2]]
  rows.append(dict(year=y,start=s,end=e,leaf_wetness_hours=r[2],mean_temperature_F=r[3],rain_inches=r[4]))
with open(out/'NEWA_2021_2025_black_rot_grouped_periods.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
def aligned(d):return d.replace(year=2000)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':16,'svg.fonttype':'none'})
fig,axes=plt.subplots(5,1,figsize=(13,10),sharex=True)
fig.subplots_adjust(left=.09,right=.855,top=.84,bottom=.16,hspace=.48)
red='#9B3C52';teal='#007F82';slate='#485467'
summary=[]
for ax,m in zip(axes,methods):
 y=int(m['Year']); rating=datetime.fromisoformat(m['Final_Assessment_Date']); rr=aligned(rating)
 ax.set_facecolor('#FFFFFF')
 first=datetime.fromisoformat(m['Experimental_Application_Dates'].split('; ')[0])
 total=sum(r['year']==y and r['end']>=first and r['start']<rating+timedelta(days=1) for r in rows)
 ax.text(1.075,.5,str(total),transform=ax.transAxes,ha='center',va='center',fontsize=25,fontweight='bold',color=slate)
 m['infection_periods_in_trial_window']=total
 ax.axvspan(rr,datetime(2000,9,5),color='#EEF0F3',zorder=0)
 for r in rows:
  if r['year']!=y:continue
  s,e=aligned(r['start']),aligned(r['end'])
  ax.plot([s,e],[.52,.52],color=red,lw=10,solid_capstyle='butt',alpha=.35 if r['start'].date()>rating.date() else 1)
 for i,d in enumerate(m['Experimental_Application_Dates'].split('; ')):
  dt=aligned(datetime.fromisoformat(d));ax.vlines(dt,.08,.72,color=teal,lw=1.1)
  ax.plot(dt,.76,marker='v',color=teal,ms=6)
  ax.text(dt,1.40 if i%2==0 else 1.03,dt.strftime('%b %-d'),ha='center',va='center',fontsize=16,color=teal)
 ax.axvline(rr,color=slate,ls=(0,(3,3)),lw=1.3)
 ax.text(rr + timedelta(days=-1.5 if y==2021 else 1.5),-.08,'Rated '+rr.strftime('%b %-d'),ha='right' if y==2021 else 'left',va='top',fontsize=16,color=slate)
 ax.set_ylim(-.45,1.6);ax.set_yticks([]);ax.set_ylabel(str(y),rotation=0,labelpad=25,fontsize=20,fontweight='bold',va='center')
 ax.spines[['left','right','top','bottom']].set_visible(False)
 ax.grid(axis='x',color='#DFE3E7',lw=.6,zorder=0);ax.tick_params(axis='x',length=0,labelbottom=False)
 def count(a,b):return sum(r['year']==y and a<=r['start']<b for r in rows)
 summary.append({'year':y,'groups_starting_May1_July14':count(datetime(y,5,1),datetime(y,7,15)),'groups_starting_June1_July14':count(datetime(y,6,1),datetime(y,7,15)),'final_rating':rating.date().isoformat(),'days_after_final_spray':m['Days_Final_Application_to_Assessment']})
axes[-1].tick_params(axis='x',labelbottom=True,pad=8)
axes[-1].set_xlim(datetime(2000,4,15),datetime(2000,9,5))
axes[-1].xaxis.set_major_locator(md.DayLocator(bymonthday=[1]));axes[-1].xaxis.set_major_formatter(md.DateFormatter('%b %d'))
fig.text(.09,.965,'Black rot infection periods, sprays, and final ratings',fontsize=23,fontweight='bold',color='#20313B')
fig.text(.925,.865,'Trial-window\nperiods',ha='center',fontsize=16,color=slate)
fig.text(.09,.928,'2021–2025  |  Geneva (AgriTech Crittenden), NY  |  Shared calendar axis',fontsize=16,color=slate)
fig.legend(handles=[Line2D([0],[0],color=red,lw=7,label='NEWA black rot period'),Line2D([0],[0],color=teal,marker='v',lw=1,label='Experimental spray'),Line2D([0],[0],color=slate,ls='--',label='Final disease rating'),Patch(facecolor='#EEF0F3',label='After final rating')],loc='lower left',bbox_to_anchor=(.085,.005),ncol=2,frameon=False,fontsize=16)
for ext in ['png','svg','pdf']:fig.savefig(out/f'NEWA_2021_2025_stacked_timeline_counts.{ext}',dpi=220,facecolor='white')
with open(out/'NEWA_common_window_summary.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=summary[0]);w.writeheader();w.writerows(summary)
print(json.dumps(summary,indent=2))
print('2025 Apr27-Aug1:',sum(r['year']==2025 and datetime(2025,4,27)<=r['start']<datetime(2025,8,2) for r in rows))

with open(out/'NEWA_trial_window_counts.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=['Year','Experimental_Application_Dates','Final_Assessment_Date','infection_periods_in_trial_window'],extrasaction='ignore');w.writeheader();w.writerows(methods)
print([(m['Year'],m['infection_periods_in_trial_window']) for m in methods])
