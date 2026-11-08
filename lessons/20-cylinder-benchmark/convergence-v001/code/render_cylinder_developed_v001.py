from pathlib import Path
import argparse,datetime,hashlib,json,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from plot_vortex_shedding_v002 import fluent_mesh,openfoam_mesh,native_field,reconstructed_vorticity
from wing_fluent_v001 import available_memory
S=Path(__file__).resolve().parent;R=Path('E:/CFD-Learning-Lab')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',required=True,type=Path);ap.add_argument('--preview',action='store_true');args=ap.parse_args();out=args.out.resolve();assert out.is_relative_to(R/'exports/figures') and not out.exists()
 assert available_memory()>=4*2**30 and shutil.disk_usage(R).free>=10*2**30
 states={};sources={};bindings={}
 for solver in ['fluent','openfoam']:
  folder=S/('analysis-'+solver+'-v001');audit=folder/('combined-'+solver+'-audit.json');a=json.loads(audit.read_text());assert a['all_new_steps_pass'] and a['steps']==800 and a['unique_native_states']==801
  for p,h in a['source_sha256'].items():assert sha(Path(p))==h,p
  sources[str(audit)]=sha(audit);m=folder/(solver+'-snapshots.json');states[solver]=json.loads(m.read_text());assert len(states[solver])==801;sources[str(m)]=sha(m);bindings[solver]=dict(audit=str(audit),audit_sha256=sha(audit),manifest=str(m),manifest_sha256=sha(m))
  for k,r in enumerate(states[solver]):assert abs(r['time_s']-(4+k*.005))<1e-8
 mesh={'fluent':fluent_mesh(Path(states['fluent'][0]['case'])),'openfoam':openfoam_mesh(Path(states['openfoam'][0]['case']))}
 chosen=[0,400,800] if args.preview else list(range(0,801,5));times=[states['fluent'][k]['time_s'] for k in chosen];out.mkdir();records=[]
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
 for frame,k in enumerate(chosen):
  t=states['fluent'][k]['time_s'];fields={solver:native_field(solver,dict(path=Path(states[solver][k]['data'])),mesh[solver]) for solver in states};vort={solver:reconstructed_vorticity(fields[solver],mesh[solver]) for solver in states}
  record=dict(time_s=t,step=round(t/.005),fields={})
  f,axes=plt.subplots(4,1,figsize=(14.4,13.2),layout='constrained')
  for j,(solver,kind) in enumerate([('fluent','speed'),('openfoam','speed'),('fluent','vorticity'),('openfoam','vorticity')]):
   ax=axes[j];name='ANSYS Fluent' if solver=='fluent' else 'OpenFOAM Foundation 14';values=fields[solver]['speed'] if kind=='speed' else vort[solver];limit=3. if kind=='speed' else 60.;assert np.isfinite(values).all()
   if kind=='speed':assert values.min()>=0 and values.max()<=limit
   im=ax.tripcolor(mesh[solver]['triangulation'],facecolors=values,cmap='viridis' if kind=='speed' else 'RdBu_r',vmin=0 if kind=='speed' else -limit,vmax=limit,rasterized=True)
   ax.set(xlim=(0,2.2),ylim=(0,.41),xlabel='x (m)',ylabel='y (m)',title=name+' | '+('native speed' if kind=='speed' else 'signed vorticity from native velocity'));ax.set_aspect('equal')
   f.colorbar(im,ax=ax,label='Speed (m/s)' if kind=='speed' else 'Vorticity (1/s)',fraction=.026,pad=.015,extend='neither' if kind=='speed' else 'both')
   record['fields'][solver+'_'+kind]=dict(minimum=float(values.min()),maximum=float(values.max()),saturated_cell_fraction=float(np.mean(np.abs(values)>limit)))
  f.suptitle('Cylinder wake | physical time t = %.3f s'%t,fontsize=18)
  f.supxlabel('21,598 cells | Re = 100 | dt = 0.005 s | fixed scales across both codes and all frames\n800 new steps (4–8 s) meet their declared numerical checks; inherited startup remains. No phase alignment.\nVorticity is an area-weighted velocity-gradient reconstruction, clipped at +/-60 1/s; native values remain unchanged.\nActual native states only. No generated intermediate motion. Benchmark matching and physical validation remain open.',fontsize=9)
  f.savefig(out/f'frame-{frame:03d}.png',dpi=100)
  if frame in {0,len(chosen)//2,len(chosen)-1}:f.savefig(out/f'state-{frame:03d}.pdf')
  plt.close(f);records.append(record)
  status=dict(state='rendering',completed_frames=frame+1,total_frames=len(chosen),time_s=t,observed_at=datetime.datetime.now().astimezone().isoformat());(out/'render-status.json').write_text(json.dumps(status,indent=2)+'\n')
 audit=dict(all_saved_field_checks_pass=True,preview=args.preview,selected_case_binding=bindings,native_source_audits_sha256=sources,records=records,plot=dict(times_s=times,speed_scale=[0,3],vorticity_scale=[-60,60],native_topology=True,physical_validation=False),rendering=dict(frames=len(chosen),frame_duration_ms=50,native_state_spacing_s=2. if args.preview else .025,no_generated_intermediate_states=True),rendered_png_sha256={str(p):sha(p) for p in sorted(out.glob('frame-*.png'))},source_sha256=sha(Path(__file__)))
 (out/'native-field-audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(dict(frames=len(chosen),out=str(out),preview=args.preview)))
if __name__=='__main__':main()