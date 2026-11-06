"""Plot measured Fluent samples and domain comparisons using Matplotlib Agg."""
from pathlib import Path
import csv,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.facecolor':'white','axes.grid':True,'grid.alpha':.18})
history=list(csv.DictReader((root/'results/lesson01-cylinder-replay/force-history.csv').open()))
x=[int(row['additional_iteration']) for row in history]
cd=[float(row['drag_coefficient']) for row in history]
cl=[float(row['lift_coefficient']) for row in history]
fig,axes=plt.subplots(1,2,figsize=(11.5,4.8),layout='constrained')
axes[0].plot(x,[(v-cd[0])/cd[0]*100 for v in cd],color='#135a7d',lw=2)
axes[0].set(xlabel='Additional solver iterations',ylabel='Change in Cd from saved value (%)',title='Drag changes very little')
axes[0].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
axes[1].plot(x,cl,color='#a14027',lw=2)
axes[1].set(xlabel='Additional solver iterations',ylabel='Lift coefficient',title='Residual asymmetry decays')
axes[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
fig.suptitle('Lesson 1 · actual Fluent force-report samples',fontsize=16,fontweight='bold')
fig.supxlabel('101 samples at report precision; iterative stability does not establish physical validation.',fontsize=10)
fig.savefig(root/'images/lesson01-force-history.png',dpi=180)
plt.close(fig)

baseline=json.loads((root/'results/lesson01-cylinder-replay/summary.json').read_text())
rows=[dict(outer_radius_m=20,drag_coefficient=baseline['final']['drag_coefficient'],cells=98304)]
for path in (root/'results').glob('lesson02-radius*/checks.json'):
    check=json.loads(path.read_text())
    assert check['numerical_checks_passed']
    geometry=json.loads((path.parent/'geometry-setup.json').read_text())
    rows.append(dict(outer_radius_m=check['outer_radius_m'],drag_coefficient=check['final_drag_coefficient'],cells=geometry['cells']))
rows.sort(key=lambda row:row['outer_radius_m'])
for i,row in enumerate(rows):
    row['change_from_preceding_percent']=None if i==0 else abs(row['drag_coefficient']-rows[i-1]['drag_coefficient'])/abs(row['drag_coefficient'])*100
source=dict(rows=rows,learning_target_percent=.2,physical_accuracy_bound=False,
            note='Same original core grid through radius 20D; outer radial rings added. No benchmark validation.')
(root/'results/lesson02-domain-comparison.json').write_text(json.dumps(source,indent=2)+'\n')
fig,axes=plt.subplots(1,2,figsize=(11.5,4.8),layout='constrained')
axes[0].plot([r['outer_radius_m'] for r in rows],[r['drag_coefficient'] for r in rows],'-o',color='#135a7d',lw=2)
axes[0].set_xscale('log',base=2)
axes[0].set_xticks([r['outer_radius_m'] for r in rows],labels=[str(r['outer_radius_m']) for r in rows])
axes[0].set(xlabel='Outer-boundary radius / cylinder diameter',ylabel='Drag coefficient',title='Moving the boundary changes drag')
labels=[f"{rows[i-1]['outer_radius_m']}→{row['outer_radius_m']}" for i,row in enumerate(rows) if i]
axes[1].bar(labels,[r['change_from_preceding_percent'] for r in rows[1:]],color='#3f7b8c',width=.6)
axes[1].axhline(.2,color='#a14027',linestyle='--',label='0.2% learning target')
axes[1].set(xlabel='Successive radius pair (D)',ylabel='Cd change / larger-domain Cd (%)',title='Test sensitivity against a declared target')
axes[1].legend(frameon=False)
fig.suptitle('Lesson 2 · cylinder domain-size sensitivity at Re = 20',fontsize=16,fontweight='bold')
fig.supxlabel('Fine inner grid preserved. A small change in one output is not an error bound or physical validation.',fontsize=10)
fig.savefig(root/'images/lesson02-domain-sensitivity.png',dpi=180)
plt.close(fig)
print(json.dumps(source,indent=2))
