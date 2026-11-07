"""Check native airfoil fields and make an equal-scale incidence comparison."""
import argparse, csv, hashlib, json, math, sys
from pathlib import Path
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm, Normalize
from matplotlib.path import Path as PolygonPath
import numpy as np
from PIL import Image

SCOPE=(-.55, 2.55, -.55, .55)
NEAR=(-.30,1.55,-.43,.43)
LABEL={"fluent":"ANSYS Fluent","openfoam":"OpenFOAM"}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
load=lambda p:json.loads(p.read_text(encoding="utf-8-sig"))
def bind(p):return (p,sha(p))

def native_fluent(run):
    case,data=(run/("airfoil-Re100-final."+s+".h5") for s in ("cas","dat"))
    with h5py.File(case,"r") as h:
        mesh=h["meshes/1"]
        xy=mesh["nodes/coords/1"][:]
        edges=mesh["faces/nodes/1/nodes"][:].reshape(-1,2).astype("int64")-1
        c0=mesh["faces/c0/1"][:].astype("int64")-1
        c1d=mesh["faces/c1/1"];c1=c1d[:].astype("int64")-1
        begin=int(c1d.attrs["minId"][0])-1
        keys=np.unique(np.r_[c0[:,None]*len(xy)+edges,
                            c1[:,None]*len(xy)+edges[begin:begin+len(c1)]])
        ids,num=np.unique(keys//len(xy),return_counts=True)
        assert np.array_equal(ids,np.arange(109600)) and np.all(num==3)
        tri=(keys%len(xy)).reshape(-1,3)
    with h5py.File(data,"r") as h:
        f=h["results/1/phase-1/cells"]
        u,v,p=[f[key+"/1"][:] for key in ("SV_U","SV_V","SV_P")]
    return xy,tri,np.column_stack((u,v)),p*2

def case_record(name, alpha, root, old, helpers, continued):
    av,top=helpers
    initial=root/f"openfoam-alpha{alpha}"
    of=initial.with_name(initial.name+"-v002") if continued and name=="openfoam" else initial
    imported=of/"input-airfoil.msh"
    targets=load(of/"declared-comparison.json")
    assert targets["angle_degrees"]==alpha and sha(imported)==targets["mesh_source_sha256"]
    source=top.read_source_mesh(imported)
    wall=source["xy"][source["edges"][source["zones"]==4]]
    assert wall.shape==(source["patch_counts"]["airfoil"],2,2) and len(wall)>=384
    if name=="fluent":
        run=old/"coarse-solve" if alpha==0 else root/"fluent-alpha5"
        assert load(run/"geometry-setup.json")==load(of/"source-setup.json")
        reports=sorted(run.glob("drag-*.txt"))+sorted(run.glob("lift-*.txt"))
        assert len(reports)==202
        inputs=[run/n for n in ("summary.json","loaded-settings.json","geometry-setup.json","declared-targets.json",
            "force-history.csv","mass-balance.txt","solver-transcript.txt","airfoil-Re100-final.cas.h5","airfoil-Re100-final.dat.h5")]+reports
        watched=dict(map(bind,inputs+[imported,of/"declared-comparison.json",of/"source-setup.json"]))
        sxy,se,cells,_=av.read_source_mesh(imported)
        values,history,_=av.check_fluent(run,targets,sxy,se,cells,retain_failed_checks=True)
        checks=dict(values["check_results"])
        skipped=[] if alpha==0 else [checks.pop("fluent.zero_angle_lift_diagnostic")]
        detail={"Cd":values["Cd"],"Cl":values["Cl"],"residuals":values["final_scaled_residuals"],
                "iteration":values["final_iteration"],"mesh_quality_min_orthogonal":values["minimum_orthogonal_quality"],
                "relative_mass_imbalance":values["relative_mass_imbalance"],"last_window":values["force_window_samples"],
                "Cd_relative_range":values["Cd_relative_range"],"Cl_absolute_range":values["Cl_absolute_range"],
                "zero_incidence_symmetry_applicable":alpha==0,"raw_zero_test_unapplied_at_5":skipped,
                "reused_saved_solution":alpha==0}
        xy,tri,u,cp=native_fluent(run)
    else:
        assert load(of/"run-status.json")["stage"]=="finished_execution"
        endpoint=4800 if continued else 3000
        t=of/str(endpoint)
        paths=[of/n for n in ("run-status.json","solver.log","mesh-check.log","system/controlDict",
                            "system/fvSchemes","system/fvSolution","constant/physicalProperties",
                            "constant/momentumTransport","0/U","0/p")]
        files=sorted((of/"postProcessing/forceCoefficients").glob("*/forceCoeffs.dat"))
        watched=dict(map(bind,paths+files+[imported,of/"declared-comparison.json",
                    t/"U",t/"p",t/"phi"]+[of/"0/polyMesh"/x for x in ("points","faces","owner","neighbour","boundary")]))
        if continued:
            status=load(of/"continuation-status.json")
            chain=load(of/"continuation-inputs.json")
            assert status["stage"]=="finished_execution" and status["from"]==3000 and status["end"]==endpoint
            assert chain["source_start_iteration"]==3000 and chain["end_iteration"]==endpoint and chain["no_reinitialization"]
            assert sha(of/"system/controlDict")==chain["revised_control_sha256"]
            for key,digest in chain["source_manifest"].items():
                assert sha(initial/key)==digest
                watched[initial/key]=digest
                if key.startswith(("processor0/3000/","processor1/3000/","3000/")):
                    assert sha(of/key)==digest
                    watched[of/key]=digest
            for key in ("continuation-status.json","continuation-inputs.json","solver-continuation.log","reconstruction-4800.log"):
                watched[of/key]=sha(of/key)
        assert "Mesh OK." in (of/"mesh-check.log").read_text()
        av.check_openfoam_settings(of,targets)
        av.check_field_conditions(av.dictionary(t/"U"),av.dictionary(t/"p"),targets)
        geometry=top.load_case_geometry(of)
        xy,tri=geometry["xy"],geometry["triangles"]
        assert len(tri)==109600
        U=top.internal_field(t/"U","vector",109600)
        assert np.max(np.abs(U[:,2]))<1e-12
        u=U[:,:2];cp=2*top.internal_field(t/"p","scalar",109600)
        _,_,_,patch_counts=av.read_source_mesh(imported)
        flux=av.boundary_flux(t/"phi",patch_counts)
        imbalance=abs(math.fsum(flux.values()))/abs(flux["inlet"])
        residual=av.residual_history([of/"solver.log"]+([of/"solver-continuation.log"] if continued else []))
        assert [row["iteration"] for row in residual]==list(range(1,endpoint+1))
        history,_=av.force_history(of/"postProcessing/forceCoefficients")
        assert history[-1]["iteration"]==endpoint
        tail=history[-100:]; assert len(tail)==100
        cd=np.array([x["Cd"] for x in tail]);cl=np.array([x["Cl"] for x in tail])
        nr={k:residual[-1][k] for k in ("Ux","Uy","p")}
        rd=float(np.ptp(cd)/abs(cd.mean()));rl=float(np.ptp(cl))
        checks={"openfoam.residuals":max(nr.values())<=1e-7,
                "openfoam.mass_balance":imbalance<=1e-6,
                "openfoam.Cd_stability":rd<=1e-4,
                "openfoam.Cl_stability":rl<=1e-5,
                "openfoam.wall_flux":abs(flux["airfoil"])/abs(flux["inlet"])<=1e-12}
        if alpha==0:checks["openfoam.zero_angle_lift_diagnostic"]=max(abs(cl))<=1e-3
        detail={"Cd":float(cd[-1]),"Cl":float(cl[-1]),"iteration":endpoint,"residuals":nr,
                "relative_mass_imbalance":imbalance,"patch_flux_m3_s":flux,"last_window":100,
                "Cd_relative_range":rd,"Cl_absolute_range":rl,"native_geometry_verified":True,
                "zero_incidence_symmetry_applicable":alpha==0}
    assert np.isfinite(u).all() and np.isfinite(cp).all()
    points=xy[tri]
    d1,d2=points[:,1]-points[:,0],points[:,2]-points[:,0]
    cross=d1[:,0]*d2[:,1]-d1[:,1]*d2[:,0]
    assert np.all(cross!=0)
    tri=tri.copy()
    tri[cross<0]=tri[cross<0][:,[0,2,1]]
    assert all(sha(p)==v for p,v in watched.items())
    checks={k:bool(v) for k,v in checks.items()}
    detail.update(checks=checks,failed_checks=[k for k,v in checks.items() if not v])
    return dict(xy=xy,tri=tri,U=u,Cp=cp,wall=wall,report=detail,inputs=watched,history=history,area=abs(cross)/2)

def visible(f, box):
    x1,x2,y1,y2=box
    pts=f["xy"][f["tri"]]
    return ((pts[:,:,0].max(axis=1)>=x1)&(pts[:,:,0].min(axis=1)<=x2)
            &(pts[:,:,1].max(axis=1)>=y1)&(pts[:,:,1].min(axis=1)<=y2))

def wall_vertices(f):
    # The original wall edges form one closed polygon. Build its ordered loop.
    ids={}
    for a,b in f["wall"]:
        a,b=tuple(a),tuple(b)
        ids.setdefault(a,[]).append(b);ids.setdefault(b,[]).append(a)
    assert len(ids)==len(f["wall"]) and all(len(v)==2 for v in ids.values())
    v0=next(iter(ids));prev=None;v=v0;out=[v0]
    for _ in range(len(ids)):
        nxt=next(b for b in ids[v] if b!=prev)
        prev,v=v,nxt
        if v==v0:break
        out.append(v)
    assert len(out)==len(ids) and v==v0
    return np.asarray(out)

def panel(ax,f,box,key,norm,cmap):
    sel=visible(f,box)
    a=ax.tripcolor(f["xy"][:,0],f["xy"][:,1],f["tri"][sel],
                   facecolors=f[key][sel],shading="flat",norm=norm,cmap=cmap,
                   edgecolors="none",rasterized=True)
    w=wall_vertices(f)
    ax.fill(w[:,0],w[:,1],color="#202633",zorder=5)
    ax.set(xlim=box[:2],ylim=box[2:],aspect="equal",xlabel="$x/c$",ylabel="$y/c$")
    return a

def plot_grid(cases, names, out, key):
    box=NEAR if key=="Cp" else SCOPE
    label=r"Pressure coefficient $C_p$" if key=="Cp" else r"Streamwise velocity $u/U_\infty$"
    v=np.concatenate([f[key][visible(f,box)] for f in cases.values()])
    norm=TwoSlopeNorm(0,vmin=v.min(),vmax=v.max()) if key=="Cp" else Normalize(v.min(),v.max())
    fig,axes=plt.subplots(len(names),2,figsize=(11,3.1*len(names)+.95),squeeze=False,layout="constrained")
    fig.suptitle(("Airfoil incidence | pressure response" if key=="Cp" else "Airfoil incidence | velocity deficit and wake")+
                 "\nNACA 0012 · Re = 100 · native saved fields",fontsize=16,fontweight="bold")
    for r,name in enumerate(names):
        for j,alpha in enumerate((0,5)):
            f=cases[name,alpha];ax=axes[r,j]
            artist=panel(ax,f,box,key,norm,"RdBu_r" if key=="Cp" else "viridis")
            ax.set_title(f"{LABEL[name]}  |  {alpha}°",loc="left",fontweight="bold")
    fig.colorbar(artist,ax=axes,label=label,shrink=.75)
    fig.savefig(out/("airfoil-pressure.png" if key=="Cp" else "airfoil-wake.png"),dpi=200)
    plt.close(fig)

def trajectories(f):
    sel=visible(f,SCOPE)
    triangles=f["tri"][sel]
    xy=f["xy"]
    tri=mtri.Triangulation(xy[:,0],xy[:,1],triangles)
    w=np.repeat(f["area"][sel],3)
    nodes=triangles.ravel()
    denom=np.bincount(nodes,weights=w,minlength=len(xy))
    denom=np.where(denom>0,denom,1.)
    interp=[mtri.LinearTriInterpolator(tri,np.bincount(nodes,
        weights=w*np.repeat(f["U"][sel,j],3),minlength=len(xy))/denom) for j in (0,1)]
    body=PolygonPath(np.vstack((wall_vertices(f),wall_vertices(f)[0])),closed=True)
    p=np.c_[np.full(12,-.48),[-.47,-.38,-.30,-.23,-.17,-.12,.12,.17,.23,.30,.38,.47]]
    states=[p.copy()]
    def velocity(q):
        return np.column_stack([a(q[:,0],q[:,1]).filled(np.nan) for a in interp])
    dt=.01
    for _ in range(400):
        ok=np.isfinite(p).all(axis=1)
        q=p[ok]
        if len(q):
            a=velocity(q);b=velocity(q+dt*a/2);c=velocity(q+dt*b/2);d=velocity(q+dt*c)
            nxt=q+dt*(a+2*b+2*c+d)/6
            good=np.isfinite(nxt).all(axis=1)&(nxt[:,0]<=SCOPE[1])&(abs(nxt[:,1])<SCOPE[3])
            # Check the full segment against the body edges, plus its endpoint.
            p0,p1=f["wall"][:,0],f["wall"][:,1]
            e=p1-p0
            delta=nxt-q
            cross=lambda a,b:a[...,0]*b[...,1]-a[...,1]*b[...,0]
            denom=cross(delta[:,None,:],e[None,:,:])
            offset=p0[None,:,:]-q[:,None,:]
            with np.errstate(invalid="ignore",divide="ignore"):
                s=cross(offset,e[None,:,:])/denom
                t=cross(offset,delta[:,None,:])/denom
            hits=np.any((s>=0)&(s<=1)&(t>=0)&(t<=1)&(denom!=0),axis=1)
            good &= ~hits & ~body.contains_points(nxt)
            nxt[~good]=np.nan
            p[ok]=nxt
        states.append(p.copy())
    return np.asarray(states)

def movie(cases,out):
    # Both 5-degree fields; common speed scale. Tracers follow saved steady U.
    names=("fluent","openfoam")
    histories={name:trajectories(cases[name,5]) for name in names}
    v=np.concatenate([cases[name,5]["ux"][visible(cases[name,5],SCOPE)] for name in names])
    norm=Normalize(v.min(),v.max())
    fig,axes=plt.subplots(2,1,figsize=(11.6,7.5),dpi=100)
    fig.subplots_adjust(left=.085,right=.88,bottom=.13,top=.82,hspace=.54)
    fig.suptitle("How the 5° airfoil turns the flow",fontsize=19,fontweight="bold",x=.085,ha="left")
    fig.text(.085,.895,"NACA 0012  |  Re = 100  |  Tracers through steady computed fields",fontsize=11)
    fig.text(.085,.044,"Travel through the saved velocity field; no time-varying CFD solution is implied.",fontsize=10)
    tails={};points={}
    for ax,name in zip(axes,names):
        f=cases[name,5]
        artist=panel(ax,f,SCOPE,"ux",norm,"viridis")
        ax.set_title(LABEL[name],loc="left",fontweight="bold")
        tails[name]=LineCollection([],colors="#fff6d8",linewidths=1.6,zorder=6,animated=True)
        ax.add_collection(tails[name])
        points[name]=ax.scatter([],[],s=16,c="white",edgecolors="#343646",linewidths=.35,zorder=7,animated=True)
    bar=fig.add_axes([.90,.21,.016,.54]);fig.colorbar(artist,cax=bar,label=r"$u/U_\infty$")
    clock=fig.text(.73,.895,"",fontsize=10,animated=True)
    fig.canvas.draw()
    bg=fig.canvas.copy_from_bbox(fig.bbox)
    frames=[]
    for j in range(56):
        i=round(j*400/55)
        fig.canvas.restore_region(bg)
        for ax,name in zip(axes,names):
            h=histories[name]
            tail=h[max(0,i-45):i+1].transpose(1,0,2)
            tails[name].set_segments([p[np.isfinite(p).all(axis=1)] for p in tail])
            points[name].set_offsets(h[i][np.isfinite(h[i]).all(axis=1)])
            ax.draw_artist(tails[name]);ax.draw_artist(points[name])
        clock.set_text(f"Tracer travel: {i*.01:.2f} s")
        fig.draw_artist(clock)
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()))
    plt.close(fig)
    frames[0].save(out/"airfoil-5deg-actual-steady-tracers.gif",save_all=True,append_images=frames[1:],duration=60,loop=0,disposal=2)
    frames[24].save(out/"airfoil-5deg-tracer-still.png")
    with Image.open(out/"airfoil-5deg-actual-steady-tracers.gif") as im:
        assert im.n_frames==56
        seen=set()
        for i in range(im.n_frames):
            im.seek(i)
            assert im.info["duration"]==60
            seen.add(hashlib.sha256(im.convert("RGB").tobytes()).hexdigest())
        assert len(seen)==56


def wall_pressure(cases,names,out):
    # Match each wall edge to the one adjacent native cell and keep actual values.
    fig,axes=plt.subplots(1,2,figsize=(10.6,4.9),sharey=True,layout="constrained")
    fig.suptitle("Pressure next to the airfoil | Re = 100",fontsize=16,fontweight="bold")
    rows=[]
    for ax,alpha in zip(axes,(0,5)):
        angle=np.deg2rad(alpha)
        for name in names:
            f=cases[name,alpha];xy=f["xy"];tri=f["tri"];n=len(xy)
            lookup={tuple(np.round(p,10)):i for i,p in enumerate(xy)}
            assert len(lookup)==n
            wall=np.array([tuple(lookup[tuple(np.round(p,10))] for p in edge) for edge in f["wall"]])
            key=lambda e:np.min(e,axis=1)*n+np.max(e,axis=1)
            pairs=tri[:,[[0,1],[1,2],[2,0]]].reshape(-1,2)
            cells=np.repeat(np.arange(len(tri)),3)
            sk=key(pairs)
            order=np.argsort(sk)
            want=key(wall)
            first=np.searchsorted(sk[order],want,"left")
            last=np.searchsorted(sk[order],want,"right")
            assert np.all(last-first==1)
            cp=f["Cp"][cells[order[first]]]
            mid=f["wall"].mean(axis=1)-[.5,0]
            x=.5+mid[:,0]*np.cos(angle)-mid[:,1]*np.sin(angle)
            y=mid[:,0]*np.sin(angle)+mid[:,1]*np.cos(angle)
            for upper in (True,False):
                idx=np.where((y>=0)==upper)[0]
                idx=idx[np.argsort(x[idx])]
                color="#1765a0" if name=="fluent" else "#bc4e2d"
                side="upper" if upper else "lower"
                ax.plot(x[idx],cp[idx],color=color,ls="-" if upper else "--",lw=1.45,
                        label=f"{LABEL[name]} · {side}")
                rows.extend((name,alpha,side,float(x[i]),float(cp[i])) for i in idx)
        ax.set(xlim=(0,1),xlabel="Chordwise position $x'/c$",title=f"{alpha}° incidence")
        ax.axhline(0,color="#526176",lw=.7)
        ax.grid(alpha=.18)
    axes[0].set_ylabel(r"Wall-adjacent cell pressure coefficient $C_p$")
    axes[0].invert_yaxis()
    axes[1].legend(loc="best",fontsize=8,frameon=False)
    fig.savefig(out/"airfoil-near-wall-cp.png",dpi=200)
    fig.savefig(out/"airfoil-near-wall-cp.svg")
    plt.close(fig)
    with (out/"near-wall-cp.csv").open("x",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["solver","angle_deg","side","chordwise_x_over_c","adjacent_cell_Cp"])
        w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("run","baseline","helpers","out"):p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--paired",action="store_true")
    p.add_argument("--continued",action="store_true")
    a=p.parse_args()
    assert a.out.drive.upper()=="E:" and not a.out.exists()
    sys.path.insert(0,str(a.helpers))
    import analyse_openfoam_airfoil_v002 as av
    import plot_openfoam_airfoil_fields_v002 as top
    names=("fluent","openfoam") if a.paired else ("fluent",)
    cases={(n,t):case_record(n,t,a.run,a.baseline,(av,top),a.continued) for n in names for t in (0,5)}
    for f in cases.values():f["ux"]=f["U"][:,0]
    a.out.mkdir(parents=True)
    plt.rcParams.update({"font.size":10,"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False})
    plot_grid(cases,names,a.out,"Cp")
    plot_grid(cases,names,a.out,"ux")
    wall_pressure(cases,names,a.out)
    if a.paired:movie(cases,a.out)
    records={f"{n}_alpha{t}":f["report"] for (n,t),f in cases.items()}
    with (a.out/"force-summary.csv").open("x",newline="",encoding="utf-8") as fp:
        w=csv.writer(fp);w.writerow(("solver","angle_deg","Cd","Cl","iterations","all_applicable_checks_passed"))
        for (n,t),f in cases.items():
            d=f["report"];w.writerow((n,t,d["Cd"],d["Cl"],d["iteration"],not d["failed_checks"]))
    manifest={str(p):sha for f in cases.values() for p,sha in f["inputs"].items()}
    assert all(sha(Path(p))==v for p,v in manifest.items())
    (a.out/"native-input-hashes-private.json").write_text(json.dumps(manifest,indent=2)+"\n")
    (a.out/"summary.json").write_text(json.dumps({"Re":100,"cells":109600,"angle_deg":[0,5],
         "cases":records,"physical_validation":False,"mesh_independence":False,
         "movie_uses_fixed_steady_fields":a.paired,"openfoam_continued_from_3000":a.continued},indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:{i:v[i] for i in ("Cd","Cl","iteration","failed_checks")} for k,v in records.items()},indent=2))
if __name__=="__main__":
    main()

