"""Plot actual native prepared fields, explicitly as inputs rather than solved flow."""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from dam_break_common_v001 import *
R=Path('E:/CFD-Learning-Lab');P=R/'local/dam-break-over-obstacle-v001';CASE=R/'native/runs/dam-break-over-obstacle-v001/qualification/openfoam-collapse-v002';OUT=R/'exports/figures/2026-11-07-dam-break-setup-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    assert not OUT.exists();proof=P/'openfoam-native-input-check-v002.json';a=json.loads(proof.read_text());g=geometry();order=g['openfoam_to_common']
    sources=[CASE/name for name in ['0/C','0/Vc','0/alpha.water','0/p_rgh']]
    for path in sources:assert sha(path)==a['native_source_sha256'][str(path)]
    alpha_raw=native_field(CASE/'0/alpha.water',1).ravel();pr_raw=native_field(CASE/'0/p_rgh',1).ravel();alpha=np.empty(NCELLS);pr=np.empty(NCELLS);alpha[order]=alpha_raw;pr[order]=pr_raw
    rho=1+alpha*999;pg=pr-rho*G*g['centres'][:,1];ij=g['cell_ij'];images=[]
    for values in [alpha,pg]:
        field=np.ma.masked_all((NY,NX));field[ij[:,1],ij[:,0]]=values;images.append(field)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    fig,ax=plt.subplots(1,2,figsize=(12,7));fig.subplots_adjust(left=.07,right=.94,bottom=.20,top=.82,wspace=.36)
    for j in range(2):
        im=ax[j].pcolormesh(g['x_edges']*1000,g['y_edges']*1000,images[j],shading='flat',cmap='Blues' if j==0 else 'viridis',vmin=0,vmax=1 if j==0 else float(pg.max()),edgecolors='#bfc6cc' if j==0 else 'none',linewidth=.12 if j==0 else 0)
        ax[j].add_patch(Rectangle((OB_LEFT*1000,0),(OB_RIGHT-OB_LEFT)*1000,OB_HEIGHT*1000,color='#8d969e',ec='#454b50',zorder=4))
        ax[j].set(aspect='equal',xlim=[0,L*1000],ylim=[0,H*1000],xlabel='x (mm)',ylabel='y (mm)',xticks=[0,146,292,438,584],yticks=[0,146,292,438,584],title='Native mesh and water fraction' if j==0 else 'Initial hydrostatic gauge-pressure input')
        fig.colorbar(im,ax=ax[j],fraction=.045,pad=.025,label='Water volume fraction' if j==0 else 'Pressure (Pa gauge)')
    fig.suptitle('Dam break over a low obstacle | prepared inputs',x=.07,ha='left',y=.975,fontsize=20,fontweight='bold')
    fig.text(.07,.918,'Actual OpenFOAM native mesh and initial fields. Collapse flow has not yet been solved.',fontsize=12)
    fig.text(.07,.115,'9,072 nonuniform quadrilateral cells | water column 146.1 x 292 mm | tank 584 x 584 mm',fontsize=11)
    fig.text(.07,.075,'Grey block: obstacle. Mesh and cell volumes checked; exact initial water area = 0.0426612 m².',fontsize=10)
    fig.text(.07,.035,'Fluent and OpenFOAM will use matching geometry and initialization. This is a setup preview, not a simulation result.',fontsize=10)
    OUT.mkdir(parents=True);fig.savefig(OUT/'dam-break-setup.png',dpi=140);fig.savefig(OUT/'dam-break-setup.pdf');plt.close(fig)
    (OUT/'input-preview-receipt.json').write_text(json.dumps(dict(input_only=True,collapse_flow_computed=False,physical_validation=False,native_input_proof_sha256=sha(proof),native_source_sha256={str(p):sha(p) for p in sources},code_sha256=sha(Path(__file__))),indent=2)+'\n');print(OUT/'dam-break-setup.png')
if __name__=='__main__':main()