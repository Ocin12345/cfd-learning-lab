from pathlib import Path
import json
R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-sloshing-time-step-v001'
a=json.loads((R/'exports/figures/2026-11-08-sloshing-dt-full-v001/selected-native-analysis.json').read_text())
for code in ['fluent','openfoam']:
 z=a['temporal_comparisons'][code];print(code,'crossings',z['base_cycle_evidence']['positive_crossing_times_s'],z['new_cycle_evidence']['positive_crossing_times_s']);print('cycles',z['base_cycle_evidence']['complete_cycles'],z['new_cycle_evidence']['complete_cycles']);print('crossing_delta_s',z['ordinal_matched_crossing_time_difference_s'])
 p=json.loads((N/(code+'-production-v001')/'production-audit-v001.json').read_text());q=json.loads((N/(code+'-qualification-v001')/'qualification-audit-v001.json').read_text());rows=q['step_records']+p['step_records'];records=q['records']+p['records'];print('maxResidual',max(max(r['residuals'].values()) for r in rows),'maxAreaDrift',max(abs(r['water_area_m2']/.01-1) for r in records),'maxAlphaExcess',max(max(0,-r['alpha_min'],r['alpha_max']-1) for r in records),'maxCo',max(r['maximum_native_interface_courant'] for r in q['courant_records']+p['courant_records']))
print('crosscode',a['small_dt_cross_code_depth_difference'])
print('waterpeak',a['global_water_dominant_peak_speed_m_s'],'allpeak',a['global_all_mixture_peak_speed_m_s'])
