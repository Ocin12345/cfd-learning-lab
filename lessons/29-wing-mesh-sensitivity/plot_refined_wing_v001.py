"""Render both saved refined-wing fields with labels tied to the native audit."""
import hashlib,json,shutil,sys
from pathlib import Path
import plot_finite_wing_native_v002 as plot

def main():
    audit=Path(sys.argv.pop(1));s=json.loads((audit/'wing-refinement-summary.json').read_text())
    for p,h in json.loads((audit/'native-source-hashes-private.json').read_text()).items():
        with Path(p).open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==h
    for key,name in [('fluent','ANSYS Fluent'),('openfoam','OpenFOAM')]:
        if not all(s['cases'][name]['checks'].values()):plot.NAMES[key]=name+' (numerical check miss)'
    plot.main()
    out=Path(sys.argv[sys.argv.index('--out')+1])
    shutil.copyfile(__file__,out/'executed-labelled-wrapper.py')
    (out/'numerical-status-labels.json').write_text(json.dumps(dict(labels=plot.NAMES,audit_summary_sha256=hashlib.sha256((audit/'wing-refinement-summary.json').read_bytes()).hexdigest()),indent=2)+'\n')
if __name__=='__main__':main()
