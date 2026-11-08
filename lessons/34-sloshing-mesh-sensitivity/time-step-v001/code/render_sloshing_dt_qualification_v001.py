from pathlib import Path
import json,hashlib
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from sloshing_dt_common_v001 import *
from analyse_nozzle_fluent_v005 import common_mapping
from audit_sloshing_dt_qualification_v003 import field
O=R/'exports/figures/2026-11-08-sloshing-dt-qualification-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def getF(case,data,mesh):
 with np.load(case/'geometry.npz') as g:m=common_mapping(case/mesh,g)
 with h5py.File(case/data) as f:
  a=np.empty(NX*NY);a[m]=cell_values(alpha_group(f));uv=[]
  for k in ['SV_U','SV_V']:
   x=np.empty(NX*NY);x[m]=cell_values(f['results/1/phase-1/cells/'+k]);uv.append(x)
 return a.reshape(NY,NX),np.column_stack(uv).reshape(NY,NX,2)
def getO(case):
 c=field(case/'0/C',3);ij=np.rint(c[:,:2]/DX-.5).astype(int);ids=ij[:,0]+NX*ij[:,1];assert len(np.unique(ids))==NX*NY
 a=np.empty(NX*NY);a[ids]=field(case/'0.02/alpha.water',1).ravel();u=np.empty((NX*NY,2));u[ids]=field(case/'0.02/U',3)[:,:2]
 return a.reshape(NY,NX),u.reshape(NY,NX,2)
def main():
 checks={}
 for code in ['fluent','openfoam']:
  p=N/(code+'-qualification-v001')/'qualification-audit-v001.json';a=json.loads(p.read_text());assert a['numerical_pass'] and a['steps']==80
  for name,h in a['source_sha256'].items():assert sha(Path(name))==h,name
  checks[code]=a
 newF=getF(N/'fluent-qualification-v001','state-0080.dat.h5','active.cas.h5');newO=getO(N/'openfoam-qualification-v001');oldF=getF(B/'fluent','tank-0040.dat.h5','tank.cas.h5');oldO=getO(B/'openfoam')
 O.mkdir(parents=True,exist_ok=False);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'text.color':'#193245'})
 edgesx=np.arange(NX+1)*DX*1000;edgesy=np.arange(NY+1)*DX*1000;limit=float(np.ceil(max(np.linalg.norm(q[1],axis=-1).max() for q in [newF,newO,oldF,oldO])*10000)/10)
 fig,axs=plt.subplots(2,2,figsize=(12.8,7.2),layout='constrained')
 for j,(name,(a,u)) in enumerate([('ANSYS Fluent',newF),('OpenFOAM',newO)]):
  for i,values in enumerate([a,np.linalg.norm(u,axis=-1)*1000]):
   im=axs[i,j].pcolormesh(edgesx,edgesy,values,cmap='Blues' if i==0 else 'viridis',vmin=0,vmax=1 if i==0 else limit,shading='flat',rasterized=True);axs[i,j].set(aspect='equal',xlabel='x (mm)',ylabel='y (mm)',title=name+(' | water fraction' if i==0 else ' | mixture speed'))
   if j==1:fig.colorbar(im,ax=axs[i,:].tolist(),label='Water volume fraction' if i==0 else 'Speed (mm/s)',shrink=.8)
 fig.suptitle('Sloshing qualification | actual t = 0.020 s | dt = 0.00025 s',fontsize=16,weight='bold');fig.supxlabel('12,800 native cells | 80 accepted steps each | a startup qualification, not the full 1.2 s comparison',fontsize=10)
 for ext in ['png','pdf']:fig.savefig(O/('actual-qualification-fields.'+ext),dpi=120)
 plt.close(fig)
 fig,axs=plt.subplots(2,2,figsize=(12.8,8),layout='constrained');x=(np.arange(NX)+.5)*DX*1000
 for name,q,color,ls in [('Fluent dt',oldF,'#16677e','--'),('Fluent dt/2',newF,'#16677e','-'),('OpenFOAM dt',oldO,'#b65c2b','--'),('OpenFOAM dt/2',newO,'#b65c2b','-')]:axs[0,0].plot(x,q[0].sum(axis=0)*DX*1000,ls,c=color,label=name)
 delta={}
 for name,q,r,color in [('Fluent',newF,oldF,'#16677e'),('OpenFOAM',newO,oldO,'#b65c2b')]:
  d=(q[0]-r[0]).sum(axis=0)*DX*1000;axs[0,1].plot(x,d,c=color,label=name);delta[name]=dict(max_depth_difference_mm=float(abs(d).max()),rms_depth_difference_mm=float(np.sqrt(np.mean(d*d))))
 for j,code in enumerate(['fluent','openfoam']):
  rows=checks[code]['step_records'];keys=list(rows[0]['residuals'])
  for k in keys:axs[1,j].semilogy([r['time_s'] for r in rows],[max(r['residuals'][k],1e-20) for r in rows],label=k)
  axs[1,j].axhline(1e-4,ls='--',c='black',lw=1);axs[1,j].set(title=('Fluent terminal scaled' if j==0 else 'OpenFOAM final-outer initial')+' residuals',xlabel='Physical time (s)',ylabel='Native residual')
 axs[0,0].set(title='Column-integrated water depth at 0.020 s',xlabel='Native fine-column centre x (mm)',ylabel='Depth (mm)');axs[0,1].set(title='Signed dt/2 minus dt depth',xlabel='Same native fine-column centre x (mm)',ylabel='Depth difference (mm)')
 for ax in axs.ravel():ax.legend(fontsize=9);ax.grid(alpha=.2)
 fig.suptitle('Same geometry and initial fields | bounded startup evidence',fontsize=16,weight='bold');fig.supxlabel('Residual norms differ. Column depth integrates volume fraction; it is not an exact traced interface or a physical-validation metric.',fontsize=9)
 for ext in ['png','pdf']:fig.savefig(O/('qualification-numerical-evidence.'+ext),dpi=120)
 plt.close(fig)
 p=[O/'actual-qualification-fields.png',O/'qualification-numerical-evidence.png']
 for f in p:
  with Image.open(f) as im:im.load()
 a=dict(actual_native_time_s=.02,new_dt_s=DT,steps_each=80,states_each=81,full_interval_complete=False,actual_image_inspection_pending=True,measurement=delta,speed_colour_range_mm_s=[0,limit],image_sha256={str(f):sha(f) for f in p},selected_audits_sha256={str(N/(code+'-qualification-v001')/'qualification-audit-v001.json'):sha(N/(code+'-qualification-v001')/'qualification-audit-v001.json') for code in checks})
 (O/'qualification-render-audit.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(delta,indent=2))
if __name__=='__main__':main()
