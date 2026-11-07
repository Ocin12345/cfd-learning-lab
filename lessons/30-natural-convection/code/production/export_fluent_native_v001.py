from pathlib import Path
import os,json,shutil,hashlib
import ansys.fluent.core as pf
from wing_fluent_v001 import available_memory
root=Path('E:/CFD-Learning-Lab'); source=root/'native/runs/natural-convection-cavity-48-v001/fluent';out=root/'exports/figures/2026-11-07-natural-convection-20s-native-v001'
assert available_memory()>4*2**30 and shutil.disk_usage(root).free>10*2**30
assert not out.exists();out.mkdir(parents=True)
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source/'qualification-final.cas.h5',source/'qualification-final.dat.h5']}
s=None
try:
 s=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
 s.settings.file.read_case_data(file_name=str(source/'qualification-final.cas.h5'))
 g=s.settings.results.graphics;g.picture.use_window_resolution=False;g.picture.x_resolution=1200;g.picture.y_resolution=1100
 g.contour.create(name='cavity-temperature');c=g.contour['cavity-temperature'];c.field='temperature';c.surfaces_list=[];c.options.filled=True;c.options.node_values=False;c.range_options.auto_range=False;c.range_options.minimum=299.6;c.range_options.maximum=300.4;c.display()
 g.views.camera.projection='orthographic';g.views.camera.position=[.025,.025,1.];g.views.camera.target=[.025,.025,0.];g.views.camera.field=[.075,.06875]
 g.picture.save_picture(file_name=str(out/'temperature-native.png'))
 (out/'export-record.json').write_text(json.dumps(dict(source=str(source),source_hashes=hashes,contour=c.get_state(),simulation_time_s=20,solver_iterations_performed=0),indent=2))
finally:
 if s:s.exit()
assert all(hashlib.sha256((source/n).read_bytes()).hexdigest()==v for n,v in hashes.items())
print('Export complete; saved simulation unchanged')
