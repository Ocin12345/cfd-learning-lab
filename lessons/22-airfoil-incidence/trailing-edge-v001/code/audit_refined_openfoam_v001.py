"""Read-only 5-degree native OpenFOAM acceptance, retaining every numerical miss."""
from pathlib import Path
import argparse,csv,json,math
import numpy as np
import analyse_openfoam_airfoil_v002 as av
import plot_openfoam_airfoil_fields_v002 as top

def outcome(residual,flux,rows):
 if len(rows)!=100:raise ValueError('Require complete original 100-sample force window')
 cd=[r['Cd'] for r in rows];cl=[r['Cl'] for r in rows]
 if not all(math.isfinite(x) for x in list(residual.values())+list(flux.values())+cd+cl):raise ValueError('Nonfinite numerical evidence')
 if flux['inlet']>=0 or flux['outlet']<=0:raise ValueError('Wrong flux signs')
 mass=abs(math.fsum(flux.values()))/abs(flux['inlet']);rd=(max(cd)-min(cd))/abs(math.fsum(cd)/100);rl=max(cl)-min(cl)
 checks={'openfoam.residuals':max(residual.values())<=1e-7,'openfoam.mass_balance':mass<=1e-6,'openfoam.Cd_stability':rd<=1e-4,'openfoam.Cl_stability':rl<=1e-5,'openfoam.wall_flux':abs(flux['airfoil'])/abs(flux['inlet'])<=1e-12}
 return dict(checks=checks,relative_mass_imbalance=mass,Cd_relative_range=rd,Cl_absolute_range=rl)

def audit(case,output,endpoint=3000,extra_logs=None):
 if output.exists():raise FileExistsError(output)
 t=json.loads((case/'declared-comparison.json').read_text());assert t['angle_degrees']==5 and t['reynolds_number']==100
 required={'relative_boundary_mass_imbalance_target':1e-6,'initial_p_U_residual_learning_target':1e-7,'cd_relative_range_target':1e-4,'cl_absolute_range_target':1e-5,'final_force_window':100}
 for k,v in required.items():assert t[k]==v,'Changed target '+k
 assert json.loads((case/'run-status.json').read_text())['stage']=='finished_execution'
 assert 'Mesh OK.' in (case/'mesh-check.log').read_text()
 for key,path in [('mesh_source_sha256',Path(t['mesh_source_path'])),('mesh_setup_sha256',case/'source-setup.json')]:assert av.sha256(path)==t[key]
 assert av.sha256(case/'input-airfoil.msh')==t['mesh_source_sha256']
 av.check_openfoam_settings(case,t)
 end=case/str(endpoint);av.check_field_conditions(av.dictionary(end/'U'),av.dictionary(end/'p'),t)
 geo=top.load_case_geometry(case);cells=t['intended_cells'];assert cells in (119580,159300)
 u=top.internal_field(end/'U','vector',cells);p=top.internal_field(end/'p','scalar',cells);assert np.isfinite(u).all() and np.isfinite(p).all() and np.max(np.abs(u[:,2]))<1e-12
 for part in ('processor0','processor1'):
  for field in ('U','p','phi'):assert (case/part/str(endpoint)/field).is_file(),'Missing native partition endpoint'
 logs=[case/'solver.log']+[case/name for name in (extra_logs or [])]
 residual=av.residual_history(logs);assert [r['iteration'] for r in residual]==list(range(1,endpoint+1))
 forces,force_files=av.force_history(case/'postProcessing/forceCoefficients');assert forces[-1]['iteration']==endpoint
 _,_,_,counts=av.read_source_mesh(case/'input-airfoil.msh');flux=av.boundary_flux(end/'phi',counts)
 lastres={k:residual[-1][k] for k in ('Ux','Uy','p')};values=outcome(lastres,flux,forces[-100:]);checks=values.pop('checks')
 paths=list(geo['inputs'])+logs+force_files+[case/n for n in ('run-status.json','mesh-check.log','source-setup.json','preparation-provenance.json','executed-prepare-script.py','system/controlDict','system/fvSolution','system/fvSchemes','constant/physicalProperties','constant/momentumTransport','0/U','0/p')]+[end/n for n in ('U','p','phi')]+[Path(__file__),Path(av.__file__),Path(top.__file__)]
 paths += [case/part/str(endpoint)/f for part in ('processor0','processor1') for f in ('U','p','phi')]
 hashes={str(p):av.sha256(p) for p in paths}
 record=dict(source_case=str(case),cells=cells,Re=100,angle_degrees=5,iteration=endpoint,Cd=forces[-1]['Cd'],Cl=forces[-1]['Cl'],residuals=lastres,residual_definition='Maximum printed initial residual within final iteration, including all pressure correction solves',patch_flux_m3_s=flux,force_window=100,zero_incidence_symmetry_applicable=False,geometry=geo['metrics'],check_results=checks,failed_checks=[k for k,v in checks.items() if not v],numerical_checks_passed=all(checks.values()),physical_validation=False,mesh_independence_established=False,source_sha256=hashes,**values)
 output.mkdir(parents=True);(output/'checks.json').write_text(json.dumps(record,indent=2)+'\n')
 for name,rows in [('forces.csv',forces),('residuals.csv',residual)]:
  with (output/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 print(json.dumps({k:record[k] for k in ('cells','Cd','Cl','residuals','relative_mass_imbalance','Cd_relative_range','Cl_absolute_range','failed_checks','numerical_checks_passed')},indent=2));return record
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--case',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--endpoint',type=int,default=3000);a.add_argument('--extra-log',action='append');p=a.parse_args();r=audit(p.case,p.out,p.endpoint,p.extra_log)
 if not r['numerical_checks_passed']:raise SystemExit(2)
