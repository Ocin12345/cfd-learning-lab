"""Audit and visualize the two-wall cavity from saved native fields only."""
import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
from matplotlib.collections import LineCollection

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("run","geometry","out","helpers"):
        ap.add_argument("--"+name,type=Path,required=True)
    a=ap.parse_args()
    assert a.out.drive.upper()=="E:" and not a.out.exists()
    sys.path.insert(0,str(a.helpers.resolve()))
    import plot_cavity_native_v003 as viz
    import analyse_cavity_fluent_v001 as fluent
    import analyse_cavity_openfoam_v003 as foam
    import prepare_openfoam_cavity_v001 as mesh
    viz.resource_guard(render_reserve=True)
    load=lambda p:json.loads(p.read_text(encoding="utf-8-sig"))
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    f,o=a.run/"fluent",a.run/"openfoam"
    last,prior=o/"2000",o/"1800"
    base=f/"cavity-iter1800.cas.h5"; end=f/"cavity-final.cas.h5"
    dat=lambda p:p.with_suffix(".h5").with_name(p.name.replace(".cas.h5",".dat.h5"))
    poly=mesh.native_poly(o)
    deps=[Path(__file__),a.geometry,
          *[Path(mod.__file__) for mod in (viz,fluent,foam,mesh,viz.native)],
          *[a.helpers/x for x in ("analyse_nozzle_fluent_v005.py",)],
          base,end,dat(base),dat(end),f/"summary.json",f/"geometry-quality.json",
          f/"final-settings.json",f/"executed-script.py",o/"input-manifest.json",
          o/"run-status.json",o/"mesh-check.log",o/"solver.log",last/"C",
          *[d/x for d in (prior,last) for x in ("U","p","phi")],
          *[o/x for x in load(o/"input-manifest.json")["prepared_hashes"]]]
    before={str(p):sha(p) for p in deps}
    assert all(sha(o/p)==h for p,h in load(o/"input-manifest.json")["prepared_hashes"].items())
    assert load(f/"summary.json")["native_final_saved"]
    assert load(o/"run-status.json")["stage"]=="finished_execution"
    g=viz.CavityGeometry(a.geometry)
    fs=load(f/"final-settings.json")
    fq=load(f/"geometry-quality.json")
    assert fq["native_numeric_gate_pass"] and fq["minimum_cell_volume_m3"]>0
    for n in ("lid","bottom","left","right"):
        w=fs["walls"][n]["momentum"]
        moving=n in ("lid","bottom")
        assert w["shear_condition"]=="No Slip"
        assert w["wall_motion"]==("Moving Wall" if moving else "Stationary Wall")
        if moving:
            assert w["relative"]=="Absolute" and w["velocity_spec"]=="Components"
            assert [float(v["value"]) for v in w["velocity_components"]]==[1.,0.]
    assert fluent.same_mesh(base,end)>0
    fn=viz.fluent_native(end,dat(end))
    fprev=g.match(viz.fluent_native(base,dat(base)))
    assert fprev["solve_evidence"]["actual_terminal_iteration"]==1800
    assert fn["solve_evidence"]["actual_terminal_iteration"]==2000
    fc0,fc1,fflux,_=fluent.native_faces(end,dat(end))
    assert np.ptp(g.areas)<1e-14
    fm=fluent.flux_metrics(fc0,fc1,fflux,np.full(6400,1/6400))
    patches=mesh.native_patches(poly)
    oq=mesh.mesh_quality((o/"mesh-check.log").read_text(encoding="utf-8"))
    og=mesh.native_geometry(poly)
    for folder in (prior,last):
        ub=foam.boundary_blocks((folder/"U").read_text(encoding="ascii"))
        pb=foam.boundary_blocks((folder/"p").read_text(encoding="ascii"))
        for name in foam.PATCHES:
            moving=name in ("lid","bottom")
            assert re.search(r"\btype\s+"+("fixedValue" if moving else "noSlip")+r"\s*;",ub[name])
            assert re.search(r"\btype\s+zeroGradient\s*;",pb[name])
            if moving or re.search(r"\bvalue\s+",ub[name]):
                assert np.max(np.abs(foam.parse_value(ub[name],80,3,"value")-([1.,0.,0.] if moving else [0.,0.,0.])))<=1e-12
        assert all(re.search(r"\btype\s+empty\s*;",b[foam.EMPTY]) for b in (ub,pb))
    om,_=foam.flux_and_divergence((last/"phi").read_text(encoding="ascii"),patches,og)
    olog=(o/"solver.log").read_text(encoding="utf-8")
    orow=foam.residual_history(olog)[-1]
    on=viz.foam_native(o,last,o/"solver.log")
    oprev=g.match({**on,"U":foam.field((prior/"U").read_text(encoding="ascii"),6400,3,[0,1,-1,0,0,0,0])[:,:2],
                  "p":foam.field((prior/"p").read_text(encoding="ascii"),6400,1,[0,2,-2,0,0,0,0])})
    fields={"fluent":g.match(fn),"openfoam":g.match(on)}
    fre=fn["solve_evidence"]["terminal_scaled_residuals"]
    ore={k:orow[k] for k in ("Ux","Uy","p")}
    report={}
    for name,field,old,res,cm,qual,peak,rms,wall in [
        ("fluent",fields["fluent"],fprev,fre,fm,fq,fm["max_dimensionless_divergence"],fm["rms_dimensionless_divergence"],fm["boundary_absolute_face_mass_flux_sum_kg_s"]),
        ("openfoam",fields["openfoam"],oprev,ore,om,oq,om["local_divergence_scaled_max_abs"],om["local_divergence_scaled_volume_weighted_RMS"],om["closed_wall_all_face_absolute_over_U_L_depth"])]:
        grid=field["gridU"]
        delta=fluent.weighted_velocity_change(old["U"],field["U"],g.areas)
        mirrored=grid[::-1]*[1,-1]
        symmetry=float(np.linalg.norm(grid-mirrored)/np.linalg.norm(grid))
        mid=g.interpolate([.5,.5],grid)
        checks={"terminal_residuals_le_1e_6":max(res.values())<=1e-6,
                "max_native_local_divergence_le_1e_5":peak<=1e-5,
                "rms_native_local_divergence_le_1e_6":rms<=1e-6,
                "normalized_wall_absolute_flux_le_1e_10":wall<=1e-10,
                "velocity_relative_change_1800_to_2000_le_1e_3":delta<=1e-3}
        report[name]={"final_iteration":2000,"terminal_residuals":res,
                      "native_mesh_quality":qual,"continuity":cm,"velocity_relative_change":delta,
                      "mirror_symmetry_relative_L2":symmetry,"interpolated_centre_U_m_s":mid.tolist(),
                      "native_vertex_set_mapping":field["mapping"],"checks":checks,
                      "failed_checks":[k for k in checks if not checks[k]],"all_declared_checks_passed":all(checks.values())}
        field["selected_iteration"]=2000
    a.out.mkdir(parents=True)
    # Same speed scale and axes; both horizontal wall velocities are drawn above/below the domain.
    decorate=viz.decorate
    def both(ax,title):
        decorate(ax,title)
        ax.plot([0,1],[0,0],color="#c15e25",lw=2.1)
        ax.annotate("",xy=(.68,.026),xytext=(.32,.026),
                    arrowprops={"arrowstyle":"->","lw":1.4,"color":"#202734"},annotation_clip=False)
    viz.decorate=both
    # The reused renderer describes fixed-field tracers; replace the title in its output.
    original_title=viz.plt.Figure.suptitle
    def title(fig,txt,*args,**kwargs):
        if txt=="Inside a lid-driven cavity":
            txt="Two moving walls, two counter-rotating cells"
        return original_title(fig,txt,*args,**kwargs)
    viz.plt.Figure.suptitle=title
    scale=(0,max(float(f["speed"].max()) for f in fields.values()))
    fig,axes=viz.plt.subplots(1,2,figsize=(10.4,6.0),layout="constrained")
    fig.suptitle("Double-lid cavity | steady laminar flow at Re = 100",fontsize=16,fontweight="bold")
    for ax,(name,field) in zip(axes,fields.items()):
        artist=viz.collection(ax,g,field["speed"],scale,"viridis")
        ax.streamplot(g.xc,g.yc,field["gridU"][:,:,0],field["gridU"][:,:,1],
                      color="#f7fbff",density=1.15,linewidth=.75,arrowsize=.85)
        both(ax,viz.NAMES[name])
    fig.colorbar(artist,ax=axes,label="Speed |U| (m/s)",shrink=.77)
    fig.savefig(a.out/"double-lid-speed-and-streamlines.png",dpi=210)
    fig.savefig(a.out/"double-lid-speed-and-streamlines.svg")
    viz.plt.close(fig)
    lines=viz.footer_lines(report)
    curves=viz.render_centrelines(g,fields,a.out,lines)
    seeds=[[x,y] for x in (.28,.52,.74) for y in (.16,.24,.36,.64,.76,.84)]
    paths={name:[g.trace(seed,field["gridU"],total=8.) for seed in seeds] for name,field in fields.items()}
    animation=viz.render_animation(g,fields,paths,a.out,lines,scale,total=8.)
    with (a.out/"centreline-velocities.csv").open("x",newline="",encoding="utf-8") as stream:
        w=csv.writer(stream)
        w.writerow(["solver","x_or_y_m","u_at_x_half_m_s","v_at_y_half_m_s"])
        for name,curve in curves.items():
            w.writerows(zip([name]*g.ny,g.yc,curve["u_at_x_half"],curve["v_at_y_half"]))
    assert before=={str(p):sha(p) for p in deps}
    difference=float(np.linalg.norm(fields["fluent"]["U"]-fields["openfoam"]["U"])/np.linalg.norm(fields["openfoam"]["U"]))
    summary={"mesh_cells":6400,"Re":100,"wall_velocity_m_s":{"top":[1,0],"bottom":[1,0],"sides":[0,0]},
             "solvers":report,"paired_velocity_relative_L2_vs_openfoam":difference,
             "fixed_native_velocity_animation":animation,"physical_validation":False,"mesh_independence":False,
             "source_bytes_unchanged":True}
    (a.out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    (a.out/"input-hashes-private.json").write_text(json.dumps(before,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"output":str(a.out),"failed":{k:v["failed_checks"] for k,v in report.items()},
                     "U_center":{k:v["interpolated_centre_U_m_s"] for k,v in report.items()},
                     "paired_relative_L2":difference,"frames":animation["frames"]},indent=2))

if __name__=="__main__":
    main()

