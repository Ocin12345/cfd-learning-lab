"""Compare native coarse/fine column depths at matching times and physical strips."""
from pathlib import Path
import argparse,csv,hashlib,json,re
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyse_nozzle_fluent_v005 import common_mapping
R=Path('E:/CFD-Learning-Lab');P=R/'local/sloshing-refinement-v001';C=R/'local/sloshing-production-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def field(path,n):
    t=path.read_text();m=re.search(r'internalField\s+uniform\s+([^;]+);',t)
    if m:return np.full(n,float(m[1]))
    m=re.search(r'internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(',t);assert m and int(m[1])==n
    end=re.search(r'\n\s*\)\s*;',t[m.end():]);assert end
    a=np.fromstring(t[m.end():m.end()+end.start()],sep=' ');assert a.size==n;return a
def chunk(group,n):
    a=np.full(n,np.nan)
    for ds in group.values():
        lo=int(np.asarray(ds.attrs['minId']).flat[0])-1;hi=int(np.asarray(ds.attrs['maxId']).flat[0]);assert np.isnan(a[lo:hi]).all();a[lo:hi]=np.asarray(ds).ravel()
    assert np.isfinite(a).all();return a
def main():
    ap=argparse.ArgumentParser();ap.add_argument('fine_figures',type=Path);ap.add_argument('out',type=Path);args=ap.parse_args();assert not args.out.exists()
    selection=C/'accepted-pair-selection-v001.json';s=json.loads(selection.read_text());plan=json.loads((P/'PLAN.json').read_text());assert sha(selection)==plan['coarse_selection_sha256']
    audit_paths=[C/'fluent-production-audit-v001.json',C/'openfoam-production-audit-v003.json'];audits=[json.loads(p.read_text()) for p in audit_paths]
    assert sha(audit_paths[0])==s['fluent_audit_sha256'] and sha(audit_paths[1])==s['openfoam_audit_sha256']
    for a in audits:
        assert a['production_checks_pass'] and a['production_complete']
        for name,h in a['native_source_sha256'].items():assert sha(Path(name))==h,name
    fine_audit=json.loads((args.fine_figures/'native-field-audit.json').read_text());assert fine_audit['all_saved_field_checks_pass']
    for name,h in fine_audit['native_source_sha256'].items():assert sha(Path(name))==h,name
    for name,h in fine_audit['numerical_audit_sha256'].items():assert sha(Path(name))==h
    fine=json.loads((args.fine_figures/'column-depth-histories.json').read_text());t=np.arange(121)*.01;assert np.allclose(fine['times_s'],t,atol=1e-12,rtol=0)
    # Fine heights use adjacent 1.25mm strips averaged into each original 2.5mm strip.
    fh=np.asarray(fine['depth_minus_50mm']);assert fh.shape==(2,121,160);fh=fh.reshape(2,121,80,2).mean(axis=3)
    coarse=np.empty((2,121,80));sources=[];f=Path(s['fluent_native']);o=Path(s['openfoam_native'])
    with np.load(f/'geometry.npz') as g:fm=common_mapping(f/'tank.cas.h5',g)
    # Actual coarse OpenFOAM block numbering is independently bound to saved cell centres.
    text=(o/'0/C').read_text();m=re.search(r'internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(',text);assert m and int(m[1])==3200
    end=re.search(r'\n\s*\)\s*;',text[m.end():]);xy=np.fromstring(text[m.end():m.end()+end.start()].replace('(',' ').replace(')',' '),sep=' ').reshape(3200,3)[:,:2]
    ij=np.rint(xy/.0025-.5).astype(int);om=ij[:,0]+80*ij[:,1];assert len(np.unique(om))==3200 and np.max(abs(xy-(ij+.5)*.0025))<1e-10
    for k in range(121):
        step=k*20;path=f/f'tank-{step:04d}.dat.h5'
        with h5py.File(path) as h:raw=chunk(h['results/1/phase-3/cells/SV_VOF'],3200)
        a=np.empty(3200);a[fm]=raw;coarse[0,k]=a.reshape(40,80).sum(axis=0)*2.5-50;sources.append(path)
        path=o/f'{step*.0005:.8g}'/'alpha.water';raw=field(path,3200);a=np.empty(3200);a[om]=raw;coarse[1,k]=a.reshape(40,80).sum(axis=0)*2.5-50;sources.append(path)
    assert np.isfinite(coarse).all() and np.max(abs(fh[:,0]-coarse[:,0]))<1e-9
    diff=fh-coarse;x=(np.arange(80)+.5)*2.5;labels=['ANSYS Fluent','OpenFOAM 14'];colors=['#176b9c','#ce652b'];out=args.out;out.mkdir(parents=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(3,2,figsize=(12,10),layout='constrained')
    for j in range(2):
        ax[0,j].plot(x,coarse[j,-1],color=colors[0],label='80 x 40 (3,200 cells)');ax[0,j].plot(x,fh[j,-1],color=colors[1],ls='--',label='160 x 80 (12,800 cells)');ax[0,j].set(title=labels[j]+' | profile at 1.20 s',xlabel='Strip centre x (mm)',ylabel='Depth minus 50 mm (mm)',xlim=[0,200]);ax[0,j].legend(fontsize=9)
        ax[1,j].plot(t,coarse[j,:,-1],color=colors[0]);ax[1,j].plot(t,fh[j,:,-1],color=colors[1],ls='--');ax[1,j].set(title='Right-wall strip: x = 197.5 to 200 mm',xlabel='Physical time (s)',ylabel='Depth minus 50 mm (mm)',xlim=[0,1.2])
        ax[2,j].plot(t,np.max(abs(diff[j]),axis=1),color='#7544a0',label='Maximum over 80 strips');ax[2,j].plot(t,np.sqrt(np.mean(diff[j]**2,axis=1)),color='#31836a',ls='--',label='RMS over 80 strips');ax[2,j].set(title='Matched-strip fine minus coarse differences',xlabel='Physical time (s)',ylabel='Absolute depth difference (mm)',xlim=[0,1.2]);ax[2,j].legend(fontsize=9)
    for row in ax:
        low=min(a.get_ylim()[0] for a in row);high=max(a.get_ylim()[1] for a in row)
        for a in row:a.set_ylim(low,high);a.grid(alpha=.18)
    fig.suptitle('Sloshing mesh sensitivity | same physics, dt = 0.5 ms',fontsize=18)
    fig.get_layout_engine().set(rect=(0,.05,1,.96));fig.text(.04,.015,'Fine columns averaged into the same 2.5 mm strips. Two grids do not establish mesh independence or physical validation.',fontsize=10)
    fig.savefig(out/'mesh-sensitivity.png',dpi=150);fig.savefig(out/'mesh-sensitivity.pdf');plt.close(fig)
    with (out/'matched-column-depths.csv').open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solver','time_s','strip_centre_x_mm','coarse_depth_minus_50mm','fine_depth_minus_50mm','fine_minus_coarse_mm'])
        for j,label in enumerate(labels):
            for k,time in enumerate(t):
                for i,xx in enumerate(x):w.writerow([label,time,xx,coarse[j,k,i],fh[j,k,i],diff[j,k,i]])
    result=dict(times_s=t.tolist(),strip_centres_x_mm=x.tolist(),same_dt_s=.0005,same_physical_conditions=True,matched_strip_width_mm=2.5,fine_initial_matches_coarse=True,mesh_independence=False,physical_validation=False,results={label:dict(maximum_depth_difference_mm=float(abs(diff[j]).max()),space_time_rms_mm=float(np.sqrt(np.mean(diff[j]**2))),final_profile_rms_mm=float(np.sqrt(np.mean(diff[j,-1]**2))),right_strip_maximum_difference_mm=float(abs(diff[j,:,-1]).max())) for j,label in enumerate(labels)},native_source_sha256={str(p):sha(p) for p in sources},fine_figure_audit_sha256=sha(args.fine_figures/'native-field-audit.json'),renderer_sha256=sha(Path(__file__)))
    (out/'comparison-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['results'],indent=2))
if __name__=='__main__':main()
