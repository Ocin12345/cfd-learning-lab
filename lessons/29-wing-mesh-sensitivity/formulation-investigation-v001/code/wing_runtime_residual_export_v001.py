import json,numpy as np

def scalar_names(session):return list(session.fields.field_data.scalar_fields.allowed_values())
def select_svars(names):
 base=['SV_CENTROID','SV_U','SV_V','SV_W','SV_P']
 if not set(base).issubset(names):raise ValueError('Missing cell-identity/value fields')
 residual=[k for k in names if k.endswith('_R') or 'RESID' in k.upper()]
 if not residual:raise ValueError('No actual retained residual variables available')
 return base+residual

def export_runtime_residuals(session,out):
 info=session.fields.solution_variable_info;v=info.get_variables_info(zone_names=['fluid'],domain_name='mixture');names=list(v.solution_variables)
 meta=dict(native_SVAR_names=names,postprocessing_scalar_names=scalar_names(session),zone_count=int(info.get_zones_info()['fluid'].count),api='PyFluent0.43publicSVARget_data;field_data.scalar_fields.allowed_values')
 (out/'runtime-residual-catalogue.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
 assert meta['zone_count']==224972
 selected=select_svars(names);arrays={}
 for name in selected:
  values=np.asarray(session.fields.solution_variable_data.get_data(variable_name=name,zone_names=['fluid'],domain_name='mixture')['fluid']);assert np.isfinite(values).all(),name;arrays[name]=values
 np.savez_compressed(out/'runtime-residual-arrays.npz',**arrays)
 meta['exported_variables']={k:dict(shape=list(v.shape),dtype=str(v.dtype)) for k,v in arrays.items()}
 (out/'runtime-residual-catalogue.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
