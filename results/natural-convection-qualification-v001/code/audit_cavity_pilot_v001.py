"""Final native-file and saved-frame audit; creates a compact numerical figure."""
from pathlib import Path
import json,re,hashlib,ast
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import analyse_render_cavity_pilot_v003 as render
from cavity_buoyancy_common_v003 import *
OUT=render.OUT;F=render.F;O=render.O

def main():
    a=json.loads((OUT/'numerical-checks.json').read_text());assert a['strict_fluent_all_steps_pass'] and len(a['fluent_all_step_checks'])==50
    assert all(all(s.values()) for s in a['heated_checks'].values())
    assert json.loads((F/'status.json').read_text())['stage']=='closed'
    settings=json.loads((F/'heated-settings.json').read_text());r=settings['residual']['equations'];assert all(s['absolute_criteria']==(1e-9 if k=='energy' else 1e-6) for k,s in r.items())
    inputpaths=[F/'qualification-final.cas.h5',F/'qualification-final.dat.h5',F/'heated-settings.json',F/'solver-transcript.txt',F/'frames.json',O/'solver.log',O/'mesh-check.log']
    control=json.loads((F/'source-hashes.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in control.items())
    frows=json.loads((F/'frames.json').read_text());fm=frows[0]['mass_per_depth_kg_m'];ofrows=[]
    for s in frows:inputpaths.append(F/s['file']);s['relative_mass_drift']=s['mass_per_depth_kg_m']/fm-1
    for i in range(51):
        q=render.ofields(i*DT);s=render.stats(q);ofrows.append(dict(time_s=i*DT,**s));inputpaths.extend(O/f'{i*DT:g}'/n for n in ['T','U','p']+(['rho'] if i else []))
    om=ofrows[0]['mass_per_depth_kg_m']
    for s in ofrows:s['relative_mass_drift']=s['mass_per_depth_kg_m']/om-1
    checks={}
    for name,rows in [('ANSYS Fluent',frows),('OpenFOAM',ofrows)]:
        checks[name]=dict(saved_states_checked=len(rows),mass=all(abs(s['relative_mass_drift'])<1e-5 for s in rows),temperature=all(COLD-1e-4<=s['min_T_K']<=s['max_T_K']<=HOT+1e-4 for s in rows),eos=all(s['eos_relative_error']<1e-6 for s in rows),max_relative_mass_drift=max(abs(s['relative_mass_drift']) for s in rows))
        assert checks[name]['mass'] and checks[name]['temperature'] and checks[name]['eos']
    poly=next(p for p in [O/'constant/polyMesh',O/'0/polyMesh'] if p.exists());inputpaths.extend(poly/n for n in ['points','faces','owner','neighbour','boundary'])
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputpaths}
    with Image.open(OUT/'heating-pilot.gif') as im:
        assert im.n_frames==6
        for i in range(im.n_frames):im.seek(i);im.load();assert im.info['duration']==500
    for p in (ROOT/'local/natural-convection-cavity-v001').glob('*.py'):ast.parse(p.read_text(encoding='utf-8-sig'))
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    times=[s['time_s'] for s in a['fluent_all_step_checks']]
    for eq in ['continuity','x-velocity','y-velocity','energy']:
        lim=1e-9 if eq=='energy' else 1e-6
        ax[0].semilogy(times,[s['residuals'][eq]/lim for s in a['fluent_all_step_checks']],label=eq)
    ax[0].axhline(1,color='black',ls='--',lw=1,label='Acceptance limit');ax[0].set(xlabel='Physical time (s)',ylabel='Final Fluent residual / declared limit',title='Every time step checked');ax[0].legend(fontsize=8);ax[0].grid(alpha=.2)
    for name,rows in [('ANSYS Fluent',frows),('OpenFOAM',ofrows)]:ax[1].plot([s['time_s'] for s in rows],[s['relative_mass_drift']*1e6 for s in rows],label=name)
    ax[1].set(xlabel='Physical time (s)',ylabel='Closed-cavity mass change (ppm)',title='Mass from native density fields');ax[1].legend();ax[1].grid(alpha=.2)
    fig.savefig(OUT/'numerical-checks.png',dpi=180);fig.savefig(OUT/'numerical-checks.pdf');plt.close(fig)
    write(OUT/'saved-state-audit.json',dict(checks=checks,fluent_steps_meeting_strict_residual_limits=50,fluent_startup_inner_limit=100,openfoam_outer_correctors=4,source_hashes=hashes,original_fluent_control_unchanged=True,animation_frames=6,animated_simulation_times_s=[0,.2,.4,.6,.8,1],physical_validation=False,full_energy_formulation_matching_pending=True))
    print(json.dumps(checks,indent=2))
if __name__=='__main__':main()
