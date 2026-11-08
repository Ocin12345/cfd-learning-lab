"""Render only accepted actual20-step thermal qualifications. No solver launch."""
import json,sys
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from jet_energy_common_v001 import *
from analyse_heated_jet_v002 import values,field,partial,common_mapping,chunk_array,counted,TDIM,NF
from jet_common_v001 import NX,NY,DX
import re
O=R/'exports/figures/2026-11-08-heated-jet-energy-qualification-v001'
def main():
 assert not O.exists();a={n:json.loads((N/(n+'-qualification-v001')/'qualification-audit-v002.json').read_text()) for n in ['fluent','openfoam']}
 for record in a.values():
  assert record['numerical_pass']
  for p,h in record['source_sha256'].items():assert sha(p)==h,p
 selection={n:dict(case=str(N/(n+'-qualification-v001')),audit_sha256=sha(N/(n+'-qualification-v001')/'qualification-audit-v002.json'),step=20,time_s=.1,full_interval=False) for n in a};dump(S/'accepted-qualification-pair-v001.json',selection);O.mkdir(parents=True)
 f=N/'fluent-qualification-v001';o=N/'openfoam-qualification-v001'
 with np.load(PARENT/'geometry.npz') as geo:m=common_mapping(f/'active.cas.h5',geo)
 FT=values(f/'heat-0020.dat.h5','SV_T');fT=np.empty_like(FT);fT[m]=FT
 c=field(o/'1250/C',NX*NY,3,[0,1,0,0,0,0,0]);m=np.rint(c[:,1]/DX-.5).astype(int)*NX+np.rint((c[:,0]+.12)/DX-.5).astype(int);assert len(np.unique(m))==NX*NY
 OT=field(o/'0.1/T',NX*NY,1,TDIM);oT=np.empty_like(OT);oT[m]=OT
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 fig,axes=plt.subplots(2,1,figsize=(11,6.5),layout='constrained',sharex=True)
 for ax,T,label in zip(axes,[fT,oT],['ANSYS Fluent','OpenFOAM']):
  im=ax.imshow(T.reshape(NY,NX),origin='lower',extent=[-6,6,0,3],aspect='equal',vmin=300,vmax=320,cmap='inferno',interpolation='nearest');ax.set(title=label,ylabel='$y/W$');ax.plot([-6,6],[0,0],color='#f27a3d',lw=2);ax.plot([-.5,.5],[3,3],color='#16b5ba',lw=2)
 axes[-1].set_xlabel('$x/W$');fig.colorbar(im,ax=axes,shrink=.85,label='Air temperature (K)');fig.suptitle('Thermal qualification at t = 0.10 s\nSame frozen flow, dt = 0.005 s; prescribed plate remains at 320 K',size=13)
 fig.savefig(O/'native-temperature-comparison.png',dpi=170);fig.savefig(O/'native-temperature-comparison.pdf');plt.close(fig)
 fig,axes=plt.subplots(3,2,figsize=(11,10),layout='constrained')
 for n,color,label in [('fluent','#1768a6','ANSYS Fluent'),('openfoam','#cb642a','OpenFOAM')]:
  r=a[n];l=r['energy_ledger'];t=np.array(l['time_s'])
  axes[0,0].plot(t,l['cumulative_net_input_J_per_m'],color=color,label=label+' boundary input');axes[0,0].plot(t,l['stored_energy_change_J_per_m'],color=color,ls='--',label=label+' storage')
  axes[0,1].plot(t,r['patch_input_W_per_m']['plate'],color=color,label=label)
  axes[1,0].semilogy(t,np.maximum(l['step_relative_error'],1e-18),color=color,label=label);axes[1,1].semilogy(t,np.maximum(l['cumulative_relative_error'],1e-18),color=color,label=label)
  ax=axes[2,0 if n=='fluent' else 1];ax.semilogy(t,[v['residual'] for v in r['step_residuals']],color=color);ax.axhline(r['thermal_residual_limit'],color='#777777',ls=':',label='Unchanged criterion');ax.set_title(label+' native terminal residual');ax.legend(frameon=False,fontsize=8)
 for ax,title,y in [(axes[0,0],'Independent discrete inputs and storage','Energy (J m$^{-1}$)'),(axes[0,1],'Prescribed-plate heat input','Heat rate (W m$^{-1}$)'),(axes[1,0],'Every-step energy mismatch','Relative mismatch'),(axes[1,1],'Cumulative prefix mismatch','Relative mismatch')]:
  ax.set(title=title,ylabel=y);ax.legend(frameon=False,fontsize=8)
 for ax in axes[1]:ax.axhline(.005,color='#777777',ls=':',label='0.5% limit')
 for ax in axes.flat:ax.set(xlabel='Physical time (s)',xlim=(0,.1));ax.grid(alpha=.15)
 fig.suptitle('Twenty actual thermal steps: qualification evidence only\nRight-endpoint backward-Euler energy accounting; residual norms differ between codes',size=12)
 fig.savefig(O/'qualification-energy-checks.png',dpi=155);fig.savefig(O/'qualification-energy-checks.pdf');plt.close(fig)
 decoded={}
 for p in O.glob('*.png'):
  with Image.open(p) as im:im.load();assert im.width>1000 and im.height>700;decoded[p.name]=dict(size=list(im.size),sha256=sha(p))
 dump(S/'qualification-export-readback-v001.json',dict(selection=selection,figures=str(O),decoded=decoded,actual_inspection_pending=True))
 print(json.dumps(dict(figures=str(O),PNG_count=len(decoded),actual_inspection_pending=True)))
if __name__=='__main__':main()
