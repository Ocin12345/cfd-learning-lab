"""PREPARE-ONLY heated-channel production-figure scaffold.

Native pilot schemas are known. Full production adapters and source binding are
intentionally unimplemented until both actual601-state audits are accepted.
No dummy fields, inferred time slices, particle animation or zero-time heat data.
"""
from dataclasses import dataclass
from pathlib import Path
import argparse,os
os.environ['TEMP']=os.environ['TMP']='E:/CFD-Learning-Lab/temp'
os.environ['MPLCONFIGDIR']='E:/CFD-Learning-Lab/temp/matplotlib-heated-channel-production'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

STATE_COUNT=601
CELL_COUNT=2304
FRAME_INDICES=np.arange(0,601,10)
FRAME_DELAY_MS=50
T_LIMITS=(300.,350.)
BOUNDARIES=('inlet','outlet','top','bottom')
NAMES={'fluent':'ANSYS Fluent','openfoam':'OpenFOAM'}
COLOURS={'fluent':'#00879c','openfoam':'#e06421'}
PLANNED_IMAGES=('heated-channel-native-temperature.png','heated-channel-outlet-and-energy.png',
    'heated-channel-boundary-power-and-balance.png','heated-channel-decoded-GIF-contact-sheet.png')
PLANNED_GIF='heated-channel-actual-temperature.gif'
FLUENT_ACTUAL_PILOT_KEYS=('temperature_K','times_s','native_step','native_dt_s','frozen_u','frozen_v','frozen_p',
    'frozen_face_mass_flux_kg_s','physical_boundary_names','native_reference_corrected_boundary_Qin_W_per_m',
    'reconstructed_boundary_Qin_W_per_m','energy_J_per_m','storage_W_per_m','outlet_bulk_temperature_step1_K')
OPENFOAM_ACTUAL_PILOT_KEYS=('time_s','solved_time_s','T_history','frozen_reference_U','frozen_reference_p',
    'frozen_reference_native_phi','boundary_names','energy_above300_J_per_m','storage_W_per_m',
    'native_Qin_W_per_m','reconstructed_Qin_W_per_m','native_Euler_relative_defect',
    'reconstructed_PDE_relative_defect','outlet_bulk_T_K','final_linear_T_residual')

@dataclass
class NativeProduction:
    """Normalized actual audited values; adapters must never pad missing data."""
    name:str
    times_s:np.ndarray
    native_step:np.ndarray
    temperature_K:np.ndarray
    frozen_U_m_s:np.ndarray
    energy_J_per_m:np.ndarray
    solved_times_s:np.ndarray
    native_Qin_W_per_m:np.ndarray
    reconstructed_Qin_W_per_m:np.ndarray
    native_Euler_relative_defect:np.ndarray
    reconstructed_PDE_relative_defect:np.ndarray
    terminal_native_thermal_residual:np.ndarray
    outlet_bulk_times_s:np.ndarray
    outlet_bulk_T_K:np.ndarray
    source_proof:dict
    numerical_audit:dict

    def validate(self):
        assert self.name in NAMES
        assert self.times_s.shape==(STATE_COUNT,) and self.native_step.shape==(STATE_COUNT,)
        assert np.array_equal(self.native_step,np.arange(601))
        assert np.allclose(self.times_s,self.native_step*.005,atol=1e-10,rtol=0)
        assert np.all(np.diff(self.times_s)>0)
        assert self.temperature_K.shape==(STATE_COUNT,CELL_COUNT) and np.isfinite(self.temperature_K).all()
        assert self.frozen_U_m_s.shape==(CELL_COUNT,2) and np.isfinite(self.frozen_U_m_s).all()
        assert self.energy_J_per_m.shape==(STATE_COUNT,) and np.isfinite(self.energy_J_per_m).all()
        assert np.array_equal(self.solved_times_s,self.times_s[1:])
        assert self.native_Qin_W_per_m.shape==self.reconstructed_Qin_W_per_m.shape==(600,4)
        assert np.isfinite(self.native_Qin_W_per_m).all() and np.isfinite(self.reconstructed_Qin_W_per_m).all()
        for value in (self.native_Euler_relative_defect,self.reconstructed_PDE_relative_defect,self.terminal_native_thermal_residual):
            assert value.shape==(600,) and np.isfinite(value).all()
        assert self.outlet_bulk_times_s.shape==self.outlet_bulk_T_K.shape
        assert np.isfinite(self.outlet_bulk_T_K).all()
        assert all(np.any(abs(self.times_s-t)<1e-10) for t in self.outlet_bulk_times_s)
        # This checks provenance structure, not numerical target outcomes.
        assert self.source_proof.get('all_native_temperatures_independently_verified') is True
        assert self.source_proof.get('actual_all_state_frozen_U_p_and_native_face_flux_verified') is True
        assert self.source_proof.get('native_boundary_power_semantics_proved') is True
        assert self.source_proof.get('initial_heat_unsolved_and_not_padded') is True
        return self


def load_fluent_production(*args,**kwargs):
    """Pending accepted full audit:601 temperature rows and600 solved power rows.

    Pilot73/73 proves firststep/nativeh300 and total-heat semantics only. The
    production adapter must bind actual native snapshots, step/time/dt, frozen
    face mass flux and energy rows to accepted output hashes; it must not expand
    the pilot's one outlet_bulk_temperature_step1_K value into a history.
    """
    raise NotImplementedError('Await accepted actual Fluent production schema; no pilot-to-production fabrication')


def load_openfoam_production(*args,**kwargs):
    """Pending full audit and native physical-boundary gradient provenance.

    Use source-bound actual native rank/physical-face values selected by the
    accepted audit. A reconstructed serialized gradT physical boundary must not
    be relabelled direct native heat flux. Respect the audit's explicit native
    versus reconstructed classification and exact native phi/upwind T evidence.
    """
    raise NotImplementedError('Await accepted actual OpenFOAM production schema and native heat provenance')


def temperature_canvas(geometry,bundles):
    """Static canvas only; caller supplies validated actual production bundles."""
    fig=plt.figure(figsize=(12.8,6.4),dpi=100)
    fig.text(.08,.963,'Heating a channel on a frozen computed flow',fontsize=20,fontweight='bold')
    fig.text(.08,.915,'Actual native temperature evolution | inlet 300 K, walls 350 K | stationary laminar flow',fontsize=9.5,color='#425e6d')
    clock=fig.text(.08,.863,'',fontsize=10,color='#304b59')
    artists={};map_axes={};width=.84;height=(12.8*width/10)/6.4
    for bottom,(name,bundle) in zip((.585,.305),bundles.items()):
        ax=fig.add_axes([.10,bottom,width,height])
        artist=PolyCollection(geometry['nodes'][geometry['quads']],array=bundle.temperature_K[0],cmap='inferno',edgecolors='none',rasterized=True)
        artist.set_clim(*T_LIMITS);ax.add_collection(artist)
        ax.set(xlim=(0,1),ylim=(0,.1),xlabel='x (m)',ylabel='y (m)');ax.set_aspect('equal');ax.set_yticks([0,.05,.1])
        ax.set_title(NAMES[name]+' | actual native temperature',loc='left',fontsize=11,pad=8,fontweight='bold')
        # Optional stationary streamline overlay: use this bundle's actual
        # frozen native U once, on the full actual cell-centre grid only.
        # Never append invented wall samples or animate these streamlines.
        artists[name]=artist;map_axes[name]=ax
    cax=fig.add_axes([.28,.198,.5,.018]);fig.colorbar(artist,cax=cax,orientation='horizontal',label='Temperature (K) | fixed common 300-350 K scale')
    fig.text(.08,.115,'Physical time is the native solver time; playback uses 61 saved frames at 50 ms each. No temporal interpolation.',fontsize=8.3,color='#425e6d')
    # Actual numerical-status footer is inserted from unmodified accepted
    # audits; expand canvas if necessary, preserving full-width10:1 maps.
    return fig,artists,clock,map_axes


def update_actual_temperature_frame(artists,clock,bundles,frame):
    state=int(FRAME_INDICES[frame]);times=[bundle.times_s[state] for bundle in bundles.values()]
    assert max(times)-min(times)<1e-10
    for name,bundle in bundles.items():artists[name].set_array(bundle.temperature_K[state])
    clock.set_text(f'Actual physical time {times[0]:.3f} s | native thermal step {state}/600 | playback 50 ms/frame')
    # Map-only hashes may repeat if actual T settles or quantizes identically.
    # Report that result; never perturb data, add particles or force motion.
    return state


def outlet_and_energy_canvas(bundles):
    fig,axes=plt.subplots(1,2,figsize=(12.8,5.6));fig.subplots_adjust(left=.09,right=.96,top=.80,bottom=.24,wspace=.28)
    fig.suptitle('Outlet warming and stored thermal energy',x=.09,y=.96,ha='left',fontsize=19,fontweight='bold')
    for name,bundle in bundles.items():
        style='--' if name=='openfoam' else '-'
        axes[0].plot(bundle.outlet_bulk_times_s,bundle.outlet_bulk_T_K,color=COLOURS[name],ls=style,label=NAMES[name])
        axes[1].plot(bundle.times_s,bundle.energy_J_per_m,color=COLOURS[name],ls=style,label=NAMES[name])
    axes[0].set(xlabel='Actual physical time (s)',ylabel='Outlet bulk\ntemperature (K)')
    axes[1].set(xlabel='Actual physical time (s)',ylabel='Energy above 300 K\nper depth (J/m)')
    for ax in axes:ax.grid(alpha=.17);ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False)
    fig.text(.09,.13,'Outlet temperature uses actual outward convective-flux weighting. Energy uses actual native volumes and depth.',fontsize=8.5,color='#425e6d')
    return fig


def power_and_defect_canvas(bundles):
    fig,axes=plt.subplots(2,2,figsize=(12.8,9.0));fig.subplots_adjust(left=.10,right=.96,top=.82,bottom=.20,wspace=.32,hspace=.45)
    fig.suptitle('Boundary power and two distinct energy diagnostics',x=.10,y=.96,ha='left',fontsize=18,fontweight='bold')
    groups=(('Heated walls',(2,3),'#c65332'),('Inlet',(0,),'#3889ac'),('Outlet',(1,),'#54795c'))
    for col,(name,bundle) in enumerate(bundles.items()):
        for title,indices,colour in groups:
            axes[0,col].plot(bundle.solved_times_s,bundle.native_Qin_W_per_m[:,indices].sum(axis=1),color=colour,lw=1.5,label=title+' native')
            axes[0,col].plot(bundle.solved_times_s,bundle.reconstructed_Qin_W_per_m[:,indices].sum(axis=1),color=colour,ls='--',lw=1.2,label=title+' reconstructed')
        axes[0,col].set(title=NAMES[name],xlabel='Solved physical time (s)',ylabel='Signed power into fluid\nper depth (W/m)')
        style='--' if name=='openfoam' else '-'
        axes[1,0].plot(bundle.solved_times_s,bundle.native_Euler_relative_defect,color=COLOURS[name],ls=style,label=NAMES[name])
        axes[1,1].plot(bundle.solved_times_s,bundle.reconstructed_PDE_relative_defect,color=COLOURS[name],ls=style,label=NAMES[name])
    for ax,target,label in ((axes[1,0],1e-5,'Native global\nEuler defect'),(axes[1,1],.005,'Reconstructed PDE\nflux defect')):
        ax.axhline(target,color='#536b78',ls=':',lw=.8,label='Declared target');ax.set_yscale('symlog',linthresh=1e-12)
        ax.set(xlabel='Solved physical time (s)',ylabel=label)
    for ax in axes.ravel():ax.grid(alpha=.17);ax.legend(frameon=False,fontsize=8,ncol=1);ax.spines[['top','right']].set_visible(False)
    fig.text(.10,.115,'Initial unsolved heat is unavailable, not zero valid power. Every physical boundary remains separate in the CSV.',fontsize=8.4,color='#425e6d')
    return fig


def export_actual_histories(*args,**kwargs):
    """Pending full schema; planned long-form source-bound exports:

    actual-native-temperatures.csv.gz: solver, native time/step, source cell,
    coordinates, actual T (1,384,704 temperatures per code; no synthetic rows).
    actual-time-metrics.csv: all601 E values, audited bulkT only where available,
    actual frozen-U/p/face-flux changes, solved flag; initial heat stays blank.
    actual-boundary-powers.csv:600 solved times x4 physical patches per code,
    native raw/report units and reference-corrected W/m, separately reconstructed
    W/m, signed difference, native diffusive/convective terms only if established.
    actual-energy-defects.csv: actual storage, boundary sums, exact declared
    normalization and separate native/reconstructed defects and target outcomes.
    actual-native-thermal-residuals.csv: native terminal residual per physical
    step with true iteration-prefix handling, not inferred steady-loop indices.
    """
    raise NotImplementedError('Bind actual accepted production history arrays first')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.parse_args()
    raise SystemExit('PREPARE-ONLY scaffold: no production adapters or render execution until both actual full audits are accepted.')

if __name__=='__main__':main()
