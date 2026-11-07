"""Qualification with Fluent universal constant verified from installed materials.h and native rho,p,T. v001 is preserved."""
from pathlib import Path
import json
ROOT=Path('E:/CFD-Learning-Lab')
L=.05; N=24; T0=300.; HOT=300.4; COLD=299.6; P0=101325.; MW=28.97
RU=8314.47247122; R=RU/MW; CP=1005.; MU=1.8e-5; PR=.71; K=MU*CP/PR; G=9.81
RHO0=P0/(R*T0); NU=MU/RHO0; ALPHA=K/(RHO0*CP); RA=G/T0*(HOT-COLD)*L**3/(NU*ALPHA)
DT=.02; CONTROL_ITERATIONS=300; HEATED_STEPS=50; INNER=100; DEPTH=.005
WALLS=['top','bottom','left','right']
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')
def specification():return dict(model='closed ideal-gas natural-convection cavity qualification',length_m=L,cells=N*N,n=N,initial_T_K=T0,hot_K=HOT,cold_K=COLD,initial_absolute_pressure_Pa=P0,molecular_weight_kg_kmol=MW,gas_constant_J_kgK=R,cp_J_kgK=CP,mu_Pas=MU,k_W_mK=K,Pr=PR,nominal_Ra=RA,gravity_m_s2=[0,-G],control_iterations=CONTROL_ITERATIONS,heated_steps=HEATED_STEPS,dt_s=DT,fluent_max_inner=INNER,openfoam_outer_correctors=4,openfoam_pressure_correctors=2,all_walls_stationary_no_slip=True,top_bottom_adiabatic_during_heating=True,ideal_gas_not_strict_boussinesq=True,qualification_only=True,physical_validation=False,full_formulation_matching_pending=True,control_max_speed_target=1e-5,control_max_temperature_deviation_K=1e-4,closed_mass_relative_drift_target=1e-5,equation_of_state_relative_error_target=1e-6)
if __name__=='__main__':print(json.dumps(specification(),indent=2))
