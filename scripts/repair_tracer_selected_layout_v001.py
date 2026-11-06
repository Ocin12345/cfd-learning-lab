"""SOURCE ONLY until root review: copy paired-v001 and repair only selected snapshot layout in fresh paired-v002."""
from pathlib import Path
import argparse,ast,copy,datetime,hashlib,importlib.util,json,os,shutil
ROOT=Path('E:/CFD-Learning-Lab');PRIVATE=ROOT/'local/prospective-tracer-v001'
SOURCE=ROOT/'exports/figures/2026-11-02-tracer-tjunction-paired-v001';OUT=ROOT/'exports/figures/2026-11-02-tracer-tjunction-paired-v002'
PLOTTER=PRIVATE/'plot_tracer_tjunction_native_v004.py';PLOTTER_SHA='60c1bb85858beedbd6f68c08818bc620a942339ed4640a76d061bf5cabc906d2'
SELECTION=PRIVATE/'visual-input-selection-v001/selection-root-approved-v001.json';SELECTION_SHA='318c9face548c3a07a7c192fc0b1e87fe56e644b7bc8250b8ce1da3695057b7c'
SOURCE_OUTPUT_MANIFEST_SHA='4a94aea29143b38f49efa22788ce3e86c0e374745a483b8e5d7b9261c035d39c'
REVIEW=PRIVATE/'visual-layout-review-v001.json';REVIEW_SHA='89bb8d97902fa43884aefb4093b2b3929f47286b560b711c738574053c715820'
RECEIPT=PRIVATE/'figure-layout-repair-v001.json'
PREFIX='tracer-tjunction';SELECTED={PREFIX+'-selected-concentration.png',PREFIX+'-selected-concentration.svg'};METADATA={'checks.json','input-manifest.json','output-manifest.json'}
EXPECTED=SELECTED|METADATA|{PREFIX+'-'+n for n in ('actual-concentration.gif','decoded-contact-sheet.png','final-concentration.png','final-concentration.svg','transport-history.png','transport-history.svg','numerical-history.png','numerical-history.svg')}|{'actual-native-scalar-cell-fields.csv.gz','record-final-scalar-fields.csv','record-scalar-step-history.csv'}

def require(v,m):
    if not v:raise ValueError(m)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def readj(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write_new(p,v):
    with Path(p).open('x',encoding='utf-8',newline='\n') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def replace_copied_metadata(p,v):
    require(p.parent==OUT and p.name in METADATA,'Only fresh copied output metadata may be replaced')
    with p.open('w',encoding='utf-8',newline='\n') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def flat_hashes(folder):
    entries=list(folder.iterdir());require({p.name for p in entries}==EXPECTED and all(p.is_file() and not p.is_symlink() for p in entries),'Exact original16 flat export artifacts required')
    return {p.name:sha(p) for p in sorted(entries)}
def verify_source_function():
    raw=PLOTTER.read_text(encoding='utf-8');own=Path(__file__).read_text(encoding='utf-8')
    original=next(n for n in ast.parse(raw).body if isinstance(n,ast.FunctionDef) and n.name=='render_snapshots');copied=next(n for n in ast.parse(own).body if isinstance(n,ast.FunctionDef) and n.name=='render_snapshots')
    original_text=ast.get_source_segment(raw,original);require(original_text.count('bottom=.11')==1 and original_text.count('[.285,.055,.43,.012]')==1,'Pinned original snapshot constants differ')
    expected=original_text.replace('bottom=.11','bottom=.15').replace('[.285,.055,.43,.012]','[.285,.075,.43,.012]')
    require(ast.get_source_segment(own,copied)==expected,'Copied snapshot function must differ only at requested two numbers');return True

def render_snapshots(g,bundles,out):
    fig,axes=plt.subplots(4,2,figsize=(12.8,12.6));fig.subplots_adjust(left=.065,right=.95,top=.90,bottom=.15,hspace=.50,wspace=.22)
    fig.text(.065,.96,'Injection, transit, arrival and early washout',fontsize=19,fontweight='bold');fig.text(.065,.928,'Actual native concentration | common 0..1 scale | frozen laminar velocity | dt = 0.025 s',fontsize=9.5,color='#425e6d')
    for row,step in enumerate((80,240,480,800)):
        for col,(name,b) in enumerate(bundles.items()):
            ax=axes[row,col];ax.add_collection(PolyCollection(g['polygons'],array=b['c'][step],cmap=cmap(),norm=Normalize(0,1,clip=False),edgecolors='none',antialiased=False));decorate(ax,f'{NAMES[name]} | actual {step*DT:.2f} s | step {step}')
    cax=fig.add_axes([.285,.075,.43,.012]);fig.colorbar(ScalarMappable(norm=Normalize(0,1,clip=False),cmap=cmap()),cax=cax,orientation='horizontal',extend='both',label='Dimensionless c | values retained without clipping')
    fig.text(.065,.018,'2D, first-order Euler/upwind; 20 s is a finite window. Complete washout and physical validation are not established.',fontsize=8,color='#425e6d');save_figure(fig,out,PREFIX+'-selected-concentration')

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-sha256',required=True);args=parser.parse_args();require(sha(__file__)==args.source_sha256,'Reviewed repair source differs')
    require(sha(PLOTTER)==PLOTTER_SHA and sha(SOURCE/'output-manifest.json')==SOURCE_OUTPUT_MANIFEST_SHA,'Exact complete reviewed original render required');require(not OUT.exists() and not RECEIPT.exists(),'Fresh versioned repair output/private receipt required');verify_source_function()
    os.environ['TEMP']=os.environ['TMP']=os.environ['TMPDIR']=str(ROOT/'temp');os.environ['PYTHONDONTWRITEBYTECODE']='1'
    spec=importlib.util.spec_from_file_location('snapshot_layout_original_plotter_v004',PLOTTER);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);before=m.resources(3.)
    # Bind local names used by the literally copied function. The original module is never patched.
    global plt,PolyCollection,cmap,Normalize,decorate,NAMES,DT,ScalarMappable,save_figure
    plt=m.plt;PolyCollection=m.PolyCollection;cmap=m.cmap;Normalize=m.Normalize;decorate=m.decorate;NAMES=m.NAMES;DT=m.DT;ScalarMappable=m.ScalarMappable;save_figure=m.save_figure
    inputs=m.Inputs();inputs.pin(dict(path=str(Path(__file__).resolve()),sha256=args.source_sha256));inputs.pin(dict(path=str(PLOTTER),sha256=PLOTTER_SHA));inputs.pin(dict(path=str(REVIEW),sha256=REVIEW_SHA))
    selected=inputs.json(dict(path=str(SELECTION),sha256=SELECTION_SHA));require(selected['root_authorizes_offline_render'] is True and selected['actual_paired_physical_steps']==800 and selected['dt_s']==.025,'Exact accepted paired field selection required')
    for ref in selected['context_inputs']:inputs.pin(ref)
    old_hashes=flat_hashes(SOURCE);original_output=readj(SOURCE/'output-manifest.json');require(original_output['manifest_self_excluded'] is True and len(original_output['files'])==15,'Original render output manifest incomplete')
    declared=set()
    for row in original_output['files']:
        p=Path(row['path']).resolve();require(p.parent==SOURCE and p.name==row['relative_path'] and p.name not in declared and old_hashes[p.name]==row['sha256'] and p.stat().st_size==row['bytes'],'Original declared artifact/path/hash differs');declared.add(p.name)
    require(declared==EXPECTED-{'output-manifest.json'},'Original output coverage differs')
    for name,pin in old_hashes.items():inputs.pin(dict(path=str(SOURCE/name),sha256=pin))
    original_inputs=readj(SOURCE/'input-manifest.json');require(original_inputs['inputs_unchanged'] is True,'Original inputs were not stable')
    for row in original_inputs['files']:
        require(row['sha256_before']==row['sha256_after'],'Original raw input changed');inputs.pin(dict(path=row['path'],sha256=row['sha256_after']))
    old_checks=readj(SOURCE/'checks.json');require(old_checks['input_hashes_stable'] is True and old_checks['input_arrays_modified'] is False and 'figure_layout_repair' not in old_checks,'Original raw check flags differ')
    g,native,_=m.load_geometry(inputs);bundles={name:m.load_history(name,selected['solvers'][name],g,native,inputs) for name in m.NAMES}
    require(m.np.array_equal(bundles['fluent']['schedule'],bundles['openfoam']['schedule']),'Accepted paired sampled schedules differ')
    OUT.mkdir()
    for name in sorted(EXPECTED):shutil.copy2(SOURCE/name,OUT/name)
    require(flat_hashes(OUT)==old_hashes,'All16 copied source artifacts must initially match byte-for-byte')
    render_snapshots(g,bundles,OUT)
    for name in EXPECTED-SELECTED:require(sha(OUT/name)==old_hashes[name],'Snapshot rendering changed another copied artifact: '+name)
    rows=inputs.finish();after=m.resources(2.5)
    new_inputs=copy.deepcopy(original_inputs);existing={Path(row['path']).resolve():row for row in new_inputs['files']}
    for row in rows:
        p=Path(row['path']).resolve()
        if p in existing:require(existing[p]['sha256_before']==existing[p]['sha256_after']==row['sha256_after'] and existing[p]['bytes']==row['bytes'],'Original input row/hash changed')
        else:new_inputs['files'].append(row);existing[p]=row
    new_checks=copy.deepcopy(old_checks);new_checks['figure_layout_repair']=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),repair_source_sha256=args.source_sha256,source_export=str(SOURCE),source_output_manifest_sha256=SOURCE_OUTPUT_MANIFEST_SHA,source_checks_sha256=old_hashes['checks.json'],regenerated_artwork=sorted(SELECTED),subplot_bottom=dict(before=.11,after=.15),colorbar_y=dict(before=.055,after=.075),footer_y_unchanged=.018,literally_copied_function_only_two_number_changes=True,accepted_fields_geometry_timepoints_colors_text_unchanged=True,all_GIF_CSV_and_other_figures_copied_byte_identically=True,original_raw_checks_and_decoded_GIF_hashes_preserved=True,original_export_unchanged=True,resources=dict(before=before,after=after),root_layout_review_pending=True,independent_visual_and_CSV_QA_not_repeated=True)
    require({k:v for k,v in new_checks.items() if k!='figure_layout_repair'}==old_checks,'Original raw checks changed')
    replace_copied_metadata(OUT/'checks.json',new_checks);replace_copied_metadata(OUT/'input-manifest.json',new_inputs)
    files=[dict(path=str(p),relative_path=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(OUT.iterdir()) if p.name!='output-manifest.json'];require(len(files)==15 and all(Path(r['path']).parent==OUT for r in files),'Every output path must point to fresh v002')
    replace_copied_metadata(OUT/'output-manifest.json',dict(manifest_self_excluded=True,files=files));new_hashes=flat_hashes(OUT)
    require(flat_hashes(SOURCE)==old_hashes,'Original export changed during repair')
    unchanged=EXPECTED-SELECTED-METADATA;require(all(new_hashes[n]==old_hashes[n] for n in unchanged),'Copied GIF/CSVs/nonselected figures changed');require(all(new_hashes[n]!=old_hashes[n] for n in SELECTED),'Requested two snapshot artifacts did not change')
    inputs.finish()
    receipt=dict(source_sha256=args.source_sha256,source_export=str(SOURCE),repaired_export=str(OUT),original_export_all16_bytes_unchanged=True,all16_original_copy_hashes_verified=True,original_raw_checks_preserved=True,decoded_GIF_hashes_preserved=True,unchanged_artifacts=sorted(unchanged),source_artifact_sha256=old_hashes,repaired_artifact_sha256=new_hashes,only_regenerated_artwork=sorted(SELECTED),metadata_rebuilt=sorted(METADATA),all_output_manifest_paths_point_to_v002=True,exact16_export_filenames=True,new_native_or_GIF_or_CSV_computation=False,accepted_fields_unchanged=True,private_receipt_outside_export=True,root_layout_review_pending=True,independent_visual_and_CSV_QA_pending=True,physical_validation=False)
    write_new(RECEIPT,receipt);print(json.dumps(dict(output=str(OUT),receipt=str(RECEIPT),receipt_sha256=sha(RECEIPT),only_regenerated_artwork=sorted(SELECTED))))
if __name__=='__main__':main()
