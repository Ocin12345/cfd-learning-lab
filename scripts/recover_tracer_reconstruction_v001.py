"""SOURCE ONLY: root-guarded reconstructPar interruption, E-backed replica recovery/export, separate owned-loop unmount. No CFD solve."""
from pathlib import Path
import argparse,datetime,hashlib,json,os,re,select,shutil,signal,stat,subprocess,time
ROOT=Path('/mnt/e/CFD-Learning-Lab');PRIVATE=ROOT/'local/prospective-tracer-v001';NATIVE=ROOT/'native/runs/2026-11-02-tracer-tjunction-v001'
ORIGINAL=NATIVE/'openfoam-production-v001';FINAL=NATIVE/'openfoam-production-reconstructed-v002';MOUNT=ROOT/'temp/tracer-reconstruction-mount-v001';IMAGE=ROOT/'temp/tracer-reconstruction-ext4-v001.img';WORK=MOUNT/'openfoam-production-working-v001';OUT=PRIVATE/'openfoam-reconstruction-recovery-v001'
PROBE=PRIVATE/'tracer-reconstruction-filesystem-probe-v001.json';PROBE_SHA='299a333ab89992d4a12b7f5d9c22797ab3192be7cdae50d3f9d6f3b01f5f6957'
ORIGINAL_AUTH=ORIGINAL/'production-provenance/root-production-authorization.json';ORIGINAL_AUTH_SHA='377fd5d23cc7c557eaa258db4d97073be86a4c61903ef0002f813a7595194f64'
WRAPPER=PRIVATE/'run-tracer-reconstruction-recovery-v001.sh';FLOAT=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'

def need(v,m):
 if not v:raise ValueError(m)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def readj(p):return json.loads(Path(p).read_text())
def win(p):return 'E:/CFD-Learning-Lab/'+str(Path(p).relative_to(ROOT))
def bind(p):return dict(path=win(p),sha256=sha(p))
def save(p,v):
 with Path(p).open('x',encoding='ascii') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def command(argv,cwd=None,log=None):
 p=subprocess.run(argv,cwd=cwd,stdout=log or subprocess.PIPE,stderr=subprocess.STDOUT if log else subprocess.PIPE,text=True,check=False)
 need(p.returncode==0,'Postprocessing/storage command failed: '+Path(argv[0]).name+' code'+str(p.returncode));return p.stdout.strip() if log is None else p.returncode

def mounted():
 need(sha(PROBE)==PROBE_SHA and sha(ORIGINAL_AUTH)==ORIGINAL_AUTH_SHA,'Exact actual probe and original production admission required');p=readj(PROBE);need(p['success'] is True and p['mounted_kept'] is True,'Actual successful mounted storage probe required')
 need(IMAGE.resolve(strict=True)==IMAGE and stat.S_ISREG(IMAGE.lstat().st_mode) and not IMAGE.is_symlink(),'Exact real regular E image required');s=IMAGE.stat();need((s.st_dev,s.st_ino)==(p['owned_new_image']['device_id'],p['owned_new_image']['inode']),'Owned E image inode changed')
 loop=p['loop_identity']['loop_device'];back=command(['losetup','--list','--noheadings','--output','BACK-FILE',loop]);need(Path(back).resolve(strict=True)==IMAGE,'Exact owned loop backing required')
 fs=json.loads(command(['findmnt','--json','--target',str(MOUNT),'--output','TARGET,SOURCE,FSTYPE']))['filesystems'];need(len(fs)==1 and fs[0]=={'target':str(MOUNT),'source':loop,'fstype':'ext4'},'Exact E-backed ext4 loop mount required');return loop

def proc(pid):
 p=Path('/proc')/str(pid)
 if not p.exists():return None
 raw=(p/'stat').read_text();tail=raw[raw.rfind(')')+2:].split()
 return dict(pid=pid,state=tail[0],ppid=int(tail[1]),start_ticks=int(tail[19]),exe=str((p/'exe').resolve()),cwd=str((p/'cwd').resolve()),argv=(p/'cmdline').read_bytes().rstrip(b'\0').decode().split('\0'))
def live(pid,start):
 try:r=proc(pid)
 except (FileNotFoundError,ProcessLookupError):return False
 return r is not None and r['start_ticks']==start and r['state']!='Z'
def target(a):
 r=proc(a['reconstruct_pid']);need(r and r['start_ticks']==a['reconstruct_start_ticks'] and r['cwd']==str(ORIGINAL) and r['exe']==a['reconstruct_exe'] and Path(r['exe']).name=='reconstructPar' and r['exe'].startswith('/opt/openfoam14/platforms/'),'Exact original reconstructPar identity required')
 need(sha(r['exe'])==a['reconstruct_exe_sha256'] and len(r['argv'])==3 and Path(r['argv'][0]).name=='reconstructPar' and r['argv'][1:]==['-time','0.15:20'],'Only exact reconstructPar0.15:20 may be interrupted')
 parent=proc(a['wrapper_linux_pid']);need(parent and r['ppid']==parent['pid'] and parent['start_ticks']==a['wrapper_start_ticks'] and Path(parent['exe']).name=='bash' and str(PRIVATE/'run-openfoam-tracer-production-v001.sh') in parent['argv'],'Exact original Bash wrapper parent required');return r

def tree(folder):
 need(folder.is_dir() and not folder.is_symlink(),'Real native tree required');rows={}
 for p in sorted(folder.rglob('*')):
  need(not p.is_symlink(),'Native symlinks unsupported')
  if p.is_file():
   relative=p.relative_to(folder).as_posix();need(stat.S_ISREG(p.stat().st_mode),'Regular native file required');rows[relative]=sha(p)
 need(len({n.casefold() for n in rows})==len(rows),'Case-folding native path collision');return rows

def copy_tree(source,destination,manifest):
 need(not destination.exists(),'Fresh versioned case copy required');destination.mkdir()
 for name,pin in manifest.items():
  p=destination/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/name,p)
 need(tree(destination)==manifest,'Full copy byte hashes differ')

def derived(name):
 if name.startswith('recovery-provenance/'):return True
 top=name.split('/')[0]
 if not re.fullmatch(FLOAT,top):return False
 step=float(top)/.025;return 6<=round(step)<=800 and abs(step-round(step))<=1e-10

def pilot_preserved(manifest):
 pins=readj(ORIGINAL_AUTH)['pilot_file_sha256'];need(len(pins)==183,'Exact accepted183-file pilot manifest required')
 for name,pin in pins.items():
  preserved='production-provenance/pilot-overridden/'+name if name in ('system/controlDict','run-status.json') else name
  need(manifest.get(preserved)==pin,'Accepted pilot bytes differ: '+name)

def admission(args):
 path=Path(args.authorization).resolve();need(PRIVATE in path.parents and sha(path)==args.authorization_sha256,'Exact root private recovery admission hash required');a=readj(path)
 expected='unmount_owned_tracer_loop_after_verified_recovery_v001' if args.mode=='unmount-owned-loop' else 'recover_tracer_openfoam_reconstruction_on_E_ext4_v001'
 need(a.get('approved') is True and a.get('approved_by')=='root' and a.get('action')==expected and a.get('source_sha256')==sha(__file__) and a.get('wrapper_sha256')==sha(WRAPPER),'Exact root worker source admission required')
 need(a.get('original_case')==str(ORIGINAL) and a.get('recovered_case')==str(FINAL) and a.get('working_replica')==str(WORK) and a.get('physical_CFD_solves_authorized')==0 and a.get('reconstruction_command')==['reconstructPar','-time','0.15:20'],'Exact zero-CFD E-only postprocessing scope required')
 need(a.get('original_production_authorization_sha256')==ORIGINAL_AUTH_SHA and a.get('probe_receipt_sha256')==PROBE_SHA,'Exact executed original/probe ancestry required');return path,a

def recover(args,a,authpath):
 for p in (OUT,WORK,FINAL):need(not os.path.lexists(p),'Fresh recovery path required; preserve existing '+p.name)
 mounted();r=target(a);need(hasattr(os,'pidfd_open') and hasattr(signal,'pidfd_send_signal'),'Identity-safe Linux pidfd signaling required');need(os.environ.get('WM_PROJECT_VERSION')=='14','Reviewed installed OpenFOAM14 environment required')
 mem={r.split(':',1)[0]:r.split(':',1)[1].strip() for r in Path('/proc/meminfo').read_text().splitlines()};need(int(mem['MemAvailable'].split()[0])*1024>=1024**3 and shutil.disk_usage(ROOT).free>=10*1024**3,'Storage-only1GiB guest RAM/E10GiB floor; CFD/audit/render floors unchanged')
 need(a.get('storage_only_RAM_floor_bytes')==1024**3 and a.get('host_available_ram_bytes',0)>=1024**3,'Actual root host storage-only admission required');status=readj(ORIGINAL/'run-status.json');need(status['stage']=='reconstructing_exactly795_new_native_states','Original worker must already be in reconstruction-only stage')
 raw=(ORIGINAL/'production-solver.log').read_text();times=[float(v) for v in re.findall(r'^Time = ('+FLOAT+r')s?\s*$',raw,re.M)];need(len(times)==795 and all(abs(t-i*.025)<=1e-12 for i,t in zip(range(6,801),times)) and re.search(r'^End\s*$',raw,re.M),'Actual original795 native scalar steps must already be computed')
 OUT.mkdir();stage='identity_guarded_reconstructPar_SIGINT';save(OUT/'interruption-target.json',dict(target=r,wrapper_linux_pid=a['wrapper_linux_pid'],wrapper_windows_pid=117552,no_solver_signal=True))
 try:
  fd=os.pidfd_open(r['pid']);again=target(a);need(all(again[k]==r[k] for k in ('pid','start_ticks','ppid','cwd','exe','argv')),'Original reconstructPar identity changed before signal')
  try:signal.pidfd_send_signal(fd,signal.SIGINT);poll=select.poll();poll.register(fd,select.POLLIN);need(bool(poll.poll(120000)),'Graceful reconstructPar SIGINT did not finish within120s; no escalation')
  finally:os.close(fd)
  deadline=time.monotonic()+120
  while live(a['wrapper_linux_pid'],a['wrapper_start_ticks']):need(time.monotonic()<deadline,'Original wrapper did not finish/sync within120s');time.sleep(1)
  oldstatus=readj(ORIGINAL/'run-status.json');oldsync=readj(ORIGINAL/'production-sync-completed.json');need(oldstatus['stage']=='failed' and oldstatus['failed_stage']=='reconstructing_exactly795_new_native_states' and oldsync['native_commands_returned'] is True and oldsync['explicit_sync_exit_code']==0 and oldsync['wrapper_exit_code']!=0,'Preserved original interruption/explicit sync proof required')
  stage='full_original_snapshot_and_verified_E_ext4_copy';before=tree(ORIGINAL);save(OUT/'original-before.json',before);pilot_preserved(before);size=sum((ORIGINAL/n).stat().st_size for n in before);rankbytes=sum((ORIGINAL/n).stat().st_size for n in before if n.startswith(('processor0/','processor1/')));need(shutil.disk_usage(MOUNT).free>=size+rankbytes+128*1024**2,'Owned image needs whole-case plus derived-output headroom')
  copy_tree(ORIGINAL,WORK,before);replica_before=tree(WORK);save(OUT/'replica-before.json',replica_before);need(before==replica_before,'Working replica differs from original')
  meta=WORK/'recovery-provenance';meta.mkdir();shutil.copyfile(__file__,meta/'executed-worker.py');shutil.copyfile(WRAPPER,meta/WRAPPER.name);shutil.copyfile(authpath,meta/'root-recovery-authorization.json');log=meta/'reconstruction.log'
  stage='reconstructing_only_saved_ranks_on_E_backed_ext4'
  with log.open('x') as f:command([a['reconstruct_exe'],'-time','0.15:20'],cwd=WORK,log=f)
  command(['sync']);after=tree(ORIGINAL);save(OUT/'original-after.json',after);need(after==before,'Original case changed after its wrapper closed');replica_after=tree(WORK);save(OUT/'replica-after.json',replica_after);need(set(before)<=set(replica_after) and all(replica_after[n]==pin or derived(n) for n,pin in before.items()) and all(n in before or derived(n) for n in replica_after),'Recovery changed inputs outside global6..800/metadata')
  stage='exporting_complete_case_to_fresh_E_native_directory';copy_tree(WORK,FINAL,replica_after);final=tree(FINAL);save(OUT/'recovered-final.json',final);need(final==replica_after,'Complete exported case byte hashes differ');command(['sync']);need(tree(ORIGINAL)==before,'Original final bytes changed')
  pilot_preserved(final)
  manifests={key:bind(OUT/name) for key,name in [('original_before','original-before.json'),('original_after','original-after.json'),('replica_before','replica-before.json'),('replica_after','replica-after.json'),('recovered_final','recovered-final.json')]}
  receipt=dict(schema='openfoam_reconstruction_recovery_v001',original_case=win(ORIGINAL),recovered_case=win(FINAL),working_replica=str(WORK),original_production_authorization=bind(ORIGINAL_AUTH),original_wrapper_windows_pid=117552,manifests=manifests,original_status=bind(ORIGINAL/'run-status.json'),original_sync=bind(ORIGINAL/'production-sync-completed.json'),original_solver_log=bind(ORIGINAL/'production-solver.log'),recovery_source_sha256={win(Path(__file__).resolve()):sha(__file__),win(WRAPPER):sha(WRAPPER)},postprocess=dict(command=['reconstructPar','-time','0.15:20'],exit_code=0,log=bind(FINAL/'recovery-provenance/reconstruction.log'),sync_exit_code=0,native_commands_returned=True),copied_rank_files_unchanged=True,pilot183_preserved=True,original_all_bytes_unchanged=True,no_CFD_replay=True,full_results_accepted=False,full_signed_scalar_budget_analysis_pending=True,owned_loop_mounted=True,dedicated_distro_close_performed=False)
  save(OUT/'execution-receipt.json',receipt);print(json.dumps(dict(receipt=win(OUT/'execution-receipt.json'),sha256=sha(OUT/'execution-receipt.json'),full_results_accepted=False)))
 except BaseException as e:save(OUT/'failure.json',dict(stage=stage,error_type=type(e).__name__,error=str(e),all_existing_files_kept=True,no_CFD_replay=True,no_signal_escalation=True));raise

def unmount(a):
 receipt=OUT/'execution-receipt.json';need(sha(receipt)==a.get('recovery_receipt_sha256') and readj(receipt)['postprocess']['native_commands_returned'] is True,'Actual successful exported recovery receipt required');loop=mounted();need(sha(OUT/'recovered-final.json')==readj(receipt)['manifests']['recovered_final']['sha256'] and FINAL.is_dir(),'Actual final manifest/export binding required');command(['sync']);command(['umount',str(MOUNT)]);need(Path(command(['losetup','--list','--noheadings','--output','BACK-FILE',loop])).resolve()==IMAGE,'Owned loop backing changed before detach');command(['losetup','--detach',loop]);save(OUT/'owned-loop-unmount.json',dict(sync_exit_code=0,owned_mount_unmounted=True,owned_loop_detached=True,loop_device=loop,image_preserved=True,dedicated_distro_close_performed=False,recovery_receipt=bind(receipt)));print('Owned E loop unmounted/detached; image preserved. Dedicated distro closure remains separate.')
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',choices=['preflight','recover','unmount-owned-loop'],required=True);p.add_argument('--authorization',required=True);p.add_argument('--authorization-sha256',required=True);args=p.parse_args();authpath,a=admission(args)
 if args.mode=='preflight':print(json.dumps(dict(read_only=True,target=target(a),owned_loop=mounted(),CFD_replay=False)));return
 if args.mode=='recover':recover(args,a,authpath)
 else:unmount(a)
if __name__=='__main__':main()
