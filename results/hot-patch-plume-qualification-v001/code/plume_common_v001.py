from pathlib import Path
import json
ROOT=Path('E:/CFD-Learning-Lab')
W=.06; H=.10; NX=48; NY=80; NC=NX*NY; DX=W/NX; DY=H/NY
HEATER_LEFT=.0225; HEATER_RIGHT=.0375
T0=300.; HOT=301.; COLD=300.; P0=101325.; MW=28.97
RU=8314.47247122; R=RU/MW; CP=1005.; MU=1.8e-5; PR=.71; K=MU*CP/PR; G=9.81
RHO0=P0/(R*T0); NU=MU/RHO0; ALPHA=K/(RHO0*CP); RA=G/T0*(HOT-COLD)*H**3/(NU*ALPHA)
DT=.01; CONTROL_ITERATIONS=600; HEATED_STEPS=100; INNER=200; DEPTH=.005
WALLS=['top','bottom','left','right','heater']
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')
def specification():return dict(model='localized floor-heater plume, initial qualification',width_m=W,height_m=H,cells=NC,nx=NX,ny=NY,heater_x_m=[HEATER_LEFT,HEATER_RIGHT],initial_T_K=T0,hot_K=HOT,cold_K=COLD,initial_absolute_pressure_Pa=P0,molecular_weight_kg_kmol=MW,gas_constant_J_kgK=R,cp_J_kgK=CP,mu_Pas=MU,k_W_mK=K,Pr=PR,nominal_Ra_based_on_height=RA,gravity_m_s2=[0,-G],control_iterations=CONTROL_ITERATIONS,heated_steps=HEATED_STEPS,dt_s=DT,fluent_max_inner=INNER,openfoam_outer_correctors=16,openfoam_pressure_correctors=2,all_walls_stationary_no_slip=True,side_walls_and_unheated_floor_adiabatic=True,ceiling_cold=True,ideal_gas_not_strict_boussinesq=True,qualification_only=True,physical_validation=False,gravity_work_source_fluent='-9.81 [m s^-2] * Density * Velocity.y',control_max_speed_target=1e-5,control_max_temperature_deviation_K=1e-4,closed_mass_relative_drift_target=1e-5,equation_of_state_relative_error_target=1e-6)
