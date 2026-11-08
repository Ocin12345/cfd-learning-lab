from pathlib import Path
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001')
s=(S/'plot_cylinder_diagnostic_v002.py').read_text().replace('2026-11-08-cylinder-convergence-diagnostic-v002','2026-11-08-cylinder-half-second-v001')
a=s.index("variants=[");b=s.index(';audits={}',a);s=s[:a]+"variants=[('Coupled','fluent-coupled-window-v001','#167d83')]"+s[b:]
s=s.replace("a['steps']==10 and a['saved_states']==11","a['steps']==100 and a['saved_states']==101")
s=s.replace("assert audits['Coupled']['numerical_window_pass'] and not audits['PISO baseline']['numerical_window_pass'] and not audits['PISO tighter pressure solve']['numerical_window_pass']","assert audits['Coupled']['numerical_window_pass'] and audits['Coupled']['all_runtime_guards_pass']")
s=s.replace('fluent-coupled-pilot-v001','fluent-coupled-window-v001').replace('4.050 s','4.550 s').replace('physical_time_s=4.05','physical_time_s=4.55').replace('Cylinder wake: a controlled convergence diagnostic','Cylinder wake: the longer convergence check').replace('same 4 s checkpoint','continues the accepted 4.05 s checkpoint').replace('all 10 new steps','all 100 new steps')
s=s.replace("O/'cylinder-diagnostic.png'","O/'cylinder-window.png'").replace("O/'cylinder-diagnostic.pdf'","O/'cylinder-window.pdf'")
compile(s,'plotwindow','exec');(S/'plot_cylinder_half_second_v001.py').write_text(s,encoding='utf-8')