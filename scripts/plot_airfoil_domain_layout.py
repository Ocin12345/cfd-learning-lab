"""Plot actual preserved airfoil-domain geometry and interface mesh; no CFD."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import shutil

import numpy as np


READER = Path(__file__).with_name('plot_openfoam_airfoil_fields.py')
if not READER.exists():
    READER = Path(__file__).with_name('executed-source-reader.py')
_spec = importlib.util.spec_from_file_location('airfoil_layout_mesh_reader', READER)
_reader = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_reader)


def require(value, message):
    if not value:
        raise ValueError(message)


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def assert_source_prefix(source, mesh):
    require(len(mesh['xy']) > len(source['xy']) and mesh['cells'] > source['cells'],
            'The extended mesh must add nodes and cells')
    require(np.array_equal(source['xy'], mesh['xy'][:len(source['xy'])]),
            'Original node IDs or exact coordinates changed')
    require(np.array_equal(source['triangles'], mesh['triangles'][:source['cells']]),
            'Original triangle cell IDs or connectivity changed')


def edge_codes(edges, nodes):
    edges = np.sort(edges, axis=1)
    return edges[:, 0]*nodes + edges[:, 1]


def load_inputs(source_path, mesh_path):
    source = _reader.read_source_mesh(source_path)
    mesh = _reader.read_source_mesh(mesh_path)
    assert_source_prefix(source, mesh)
    old_setup_path, setup_path = source_path.parent/'setup.json', mesh_path.parent/'setup.json'
    provenance_path = mesh_path.parent/'domain-extension-provenance.json'
    old_setup, setup, provenance = [json.loads(p.read_text()) for p in
                                   [old_setup_path, setup_path, provenance_path]]
    require(sha256(source_path) == provenance['source_sha256']
            and sha256(mesh_path) == provenance['mesh_sha256']
            and sha256(old_setup_path) == provenance['source_setup_sha256'],
            'Domain-extension source/mesh/setup hashes differ')
    snapshots = [(mesh_path.parent/'executed-mesher.py', 'script_sha256'),
                 (mesh_path.parent/'executed-source-reader.py', 'source_reader_sha256')]
    for path, key in snapshots:
        require(sha256(path) == provenance[key], 'Domain preparation snapshot hash differs')
    require(source['cells'] == old_setup['cells'] and len(source['xy']) == old_setup['nodes']
            and mesh['cells'] == setup['cells'] == provenance['cells']
            and len(mesh['xy']) == setup['nodes'] == provenance['nodes'], 'Mesh counts differ')
    require(setup['domain_extension'] == provenance, 'Setup/provenance receipts differ')
    require(provenance['mesher']['method'] == 'deterministic_radial_triangles'
            and provenance['mesher']['radial_layers'] == 96
            and provenance['mesher']['gmsh_used'] is False,
            'This drawing/caption requires the actual radial-96 extension')
    chord = float(setup['chord_m'])
    require(chord > 0 and chord == old_setup['chord_m'], 'Invalid/changed reference chord')
    wall_old = source['edges'][source['zones'] == 4]
    wall_new = mesh['edges'][mesh['zones'] == 4]
    require(np.array_equal(np.sort(edge_codes(wall_old, len(mesh['xy']))),
                           np.sort(edge_codes(wall_new, len(mesh['xy'])))), 'Wall geometry changed')
    shared = source['edges'][np.isin(source['zones'], [5, 6])]
    old_ring = source['xy'][np.unique(shared)]
    outer_ring = mesh['xy'][np.unique(mesh['edges'][np.isin(mesh['zones'], [5, 6])])]
    scaled_ring = (old_ring-[.5,0])*2+[.5,0]
    distance = np.linalg.norm(outer_ring[:,None,:]-scaled_ring[None,:,:], axis=2)
    nearest = distance.argmin(axis=1)
    require(len(outer_ring) == len(old_ring) and len(np.unique(nearest)) == len(old_ring)
            and np.max(distance[np.arange(len(outer_ring)),nearest]) <= 1e-10,
            'Outer polygon/subdivision coordinates are not the exact twofold radial scale')
    shared_codes = edge_codes(shared, len(mesh['xy']))
    is_shared = np.isin(edge_codes(mesh['edges'], len(mesh['xy'])), shared_codes)
    require(int(is_shared.sum()) == len(shared) and np.all(mesh['zones'][is_shared] == 3)
            and np.all((mesh['owner'][is_shared] < source['cells'])
                       ^ (mesh['neighbour'][is_shared] < source['cells'])),
            'The retained farfield must join original and added cells as an interior interface')
    for zone, inlet in [(5, True), (6, False)]:
        edges = mesh['edges'][mesh['zones'] == zone]
        require(np.all((mesh['xy'][edges, 0].mean(axis=1) < .5) == inlet),
                'New outer patch split differs from actual x-midpoint rule')
    files = [source_path, mesh_path, old_setup_path, setup_path, provenance_path, READER,
             *(path for path, _ in snapshots)]
    checks = dict(original_node_ids_and_coordinates_exact=True, original_cell_ids_and_triangles_exact=True,
                  original_wall_edge_identity_exact=True, former_farfield_edges_now_interior=len(shared),
                  new_outer_patch_assignment_verified=True, source_and_preparation_hashes_verified=True)
    return source, mesh, old_setup, setup, provenance, is_shared, files, checks


def visible_segments(xy, edges, limits):
    segments = xy[edges]
    xmin, xmax, ymin, ymax = limits
    visible = ((segments[:, :, 0].min(axis=1) <= xmax) & (segments[:, :, 0].max(axis=1) >= xmin)
               & (segments[:, :, 1].min(axis=1) <= ymax) & (segments[:, :, 1].max(axis=1) >= ymin))
    return segments[visible]


def ordered_polyline(xy, edges):
    """Join exactly the supplied polygon edges without resetting dash phase."""
    adjacency = {}
    for a, b in edges:
        adjacency.setdefault(int(a), []).append(int(b))
        adjacency.setdefault(int(b), []).append(int(a))
    require(all(len(v) == 2 for v in adjacency.values()), 'Polygon is not a closed edge loop')
    first = min(adjacency)
    sequence, previous, current = [first], None, first
    while True:
        candidate = next(v for v in adjacency[current] if v != previous)
        sequence.append(candidate)
        if candidate == first:
            break
        require(candidate not in sequence[:-1], 'Polygon repeats a node')
        previous, current = current, candidate
    require(len(sequence) == len(edges)+1, 'Polygon has disconnected edge loops')
    return xy[sequence]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True, help='Original R20 Fluent ASCII mesh')
    parser.add_argument('--mesh', type=Path, required=True, help='Extended R40 Fluent ASCII mesh')
    parser.add_argument('--output-dir', type=Path, required=True, help='Fresh external figure directory')
    args = parser.parse_args()
    source_path, mesh_path, out = [p.resolve() for p in [args.source, args.mesh, args.output_dir]]
    root, repo = Path('E:/CFD-Learning-Lab').resolve(), Path(__file__).resolve().parents[1]
    require(all(p.is_relative_to(root) for p in [source_path, mesh_path, out]),
            'Keep source/output under the project root')
    require(not out.exists() and not out.is_relative_to(repo), 'Use a fresh external figure directory')
    source, mesh, old_setup, setup, provenance, is_shared, files, checks = load_inputs(source_path, mesh_path)
    chord = float(setup['chord_m'])
    old_radius, radius = old_setup['outer_radius_m']/chord, setup['outer_radius_m']/chord
    require(old_radius == 20 and radius == 40 and chord == 1,
            'This layout/crop is specifically for the unit-chord R20-to-R40 comparison')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    plt.rcParams.update({'font.family':'STIXGeneral', 'mathtext.fontset':'stix', 'font.size':9,
                         'axes.labelsize':10, 'axes.linewidth':.65, 'xtick.direction':'in',
                         'ytick.direction':'in', 'xtick.top':True, 'ytick.right':True,
                         'legend.fontsize':7.5, 'pdf.fonttype':42, 'svg.fonttype':'path',
                         'savefig.bbox':None})
    fig, axes = plt.subplots(1, 2, figsize=(7.05, 3.3), layout='constrained')
    fig.suptitle('Geometry and mesh only', fontsize=10.5)
    blue, orange, gray = '#0072B2', '#D55E00', '#91979C'
    xy, old_xy = mesh['xy']/chord, source['xy']/chord
    ax = axes[0]
    for zone, color in [(5, blue), (6, orange)]:
        ax.add_collection(LineCollection(xy[mesh['edges'][mesh['zones'] == zone]],
                                        colors=color, linewidths=1.2))
    old_far = source['edges'][np.isin(source['zones'], [5, 6])]
    old_polygon = ordered_polyline(old_xy, old_far)
    ax.plot(old_polygon[:,0], old_polygon[:,1], color=gray, linewidth=.9, linestyle='--')
    ax.add_collection(LineCollection(old_xy[source['edges'][source['zones'] == 4]],
                                    colors='black', linewidths=.8))
    ax.set_xlim(-43, 44); ax.set_ylim(-43, 43)
    ax.set_xticks([-40, -20, 0, 20, 40]); ax.set_yticks([-40, -20, 0, 20, 40])
    ax.annotate(r'$U_\infty$', xy=(-13, 28), xytext=(-34, 28), va='center', ha='left',
                arrowprops=dict(arrowstyle='->', color='black', lw=.8))
    ax.annotate(r'$R/c=20$', xy=(.5, 20), xytext=(7, 18), ha='left', fontsize=9,
                arrowprops=dict(arrowstyle='-', color=gray, lw=.7))
    ax.annotate(r'$R/c=40$', xy=(.5, 40), xytext=(9, 35), ha='left', fontsize=9,
                arrowprops=dict(arrowstyle='-', color=orange, lw=.7))
    ax.annotate('Airfoil', xy=(.5, 0), xytext=(-24, -10), ha='left', fontsize=8,
                arrowprops=dict(arrowstyle='-', color='black', lw=.65))
    ax.legend([Line2D([], [], color=blue, lw=1.2), Line2D([], [], color=orange, lw=1.2),
               Line2D([], [], color=gray, lw=.9, ls='--')],
              ['R40 inlet', 'R40 outlet', 'Retained R20 polygon'],
              loc='lower left', frameon=True, facecolor='white', edgecolor='none',
              framealpha=.96, handlelength=1.8, labelspacing=.25)
    zoom = (20.15, 20.85, -.35, .35)
    ax.add_patch(Rectangle((zoom[0], zoom[2]), zoom[1]-zoom[0], zoom[3]-zoom[2],
                           fill=False, edgecolor='black', linewidth=.7))
    ax.text(23, -2, '(b)', fontsize=8.5)
    ax.set_title('(a) Actual boundary polygons', loc='left', fontsize=10)
    ax = axes[1]
    old_shared = np.isin(edge_codes(source['edges'], len(mesh['xy'])), edge_codes(old_far, len(mesh['xy'])))
    retained = visible_segments(old_xy, source['edges'][~old_shared], zoom)
    added = ((mesh['owner'] >= source['cells']) | (mesh['neighbour'] >= source['cells'])) & ~is_shared
    annulus = visible_segments(xy, mesh['edges'][added], zoom)
    shared = visible_segments(old_xy, old_far, zoom)
    require(len(retained) > 0 and len(annulus) > 0 and len(shared) > 0, 'Empty interface zoom region')
    for segments, color, width, style in [(retained, gray, .65, '-'), (annulus, blue, .65, '-'),
                                         (shared, 'black', 1.0, '--')]:
        ax.add_collection(LineCollection(segments, colors=color, linewidths=width, linestyles=style))
    ax.set_xlim(*zoom[:2]); ax.set_ylim(*zoom[2:])
    ax.set_xticks([20.2, 20.4, 20.6, 20.8]); ax.set_yticks([-.3, 0, .3])
    ax.set_title('(b) Retained interface and triangles', loc='left', fontsize=10)
    ax.legend([Line2D([], [], color=gray, lw=.8), Line2D([], [], color=blue, lw=.8),
               Line2D([], [], color='black', lw=1, ls='--')],
              ['Retained mesh', 'Added annulus', 'Shared interface'], loc='upper left',
              frameon=True, facecolor='white', edgecolor='none', framealpha=.96,
              handlelength=1.8, labelspacing=.25)
    for ax in axes:
        ax.set_aspect('equal', adjustable='box'); ax.set_xlabel(r'$x/c$'); ax.set_ylabel(r'$y/c$')
    fig.canvas.draw(); fig.set_layout_engine('none')
    # Do not use tight-bbox exports: fixed-canvas PDF raster/bbox behavior differs.
    out.mkdir(parents=True)
    base = out/'fig10-airfoil-domain-layout'
    for suffix in ['png', 'pdf', 'svg']:
        fig.savefig(base.with_suffix('.'+suffix), dpi=600)
    plt.close(fig)
    caption = (
        'Figure 10. Geometry and mesh only; no flow solution is shown. '
        '(a) Actual R40 polygon edges are coloured by the Fluent velocity-inlet (blue) and '
        'pressure-outlet (orange) patches; the original R20 polygon is dashed gray. '
        'Radii refer to the original polygon vertices about (x/c,y/c)=(0.5,0); chord-subdivision '
        'points lie on straight polygon segments. The tiny black airfoil uses the actual retained '
        'wall edges, and the arrow indicates the specified positive-x inlet direction. '
        '(b) Actual triangle edges at x/c=20.15..20.85, y/c=-0.35..0.35 show retained edges in gray, '
        'added annulus edges in blue and the former outer boundary as a dashed black conformal '
        'interior interface. No triangles were interpolated, invented or remeshed for this drawing. '
        f"The R40 mesh contains {mesh['cells']:,} cells and {len(mesh['xy']):,} nodes, retaining all "
        f"{source['cells']:,} original cell IDs/triangles and {len(source['xy']):,} original node IDs/coordinates exactly. "
        'The annulus has 96 deterministic geometric radial strips and alternating triangle diagonals; '
        'its outer polygon is the exact uniform twofold radial scale of the original polygon. '
        'This figure illustrates a controlled domain extension, not a flow, convergence or physical-validation result.\n')
    (out/'CAPTION.md').write_text(caption, encoding='utf-8')
    shutil.copy2(Path(__file__), out/'executed-plotter.py')
    shutil.copy2(READER, out/'executed-source-reader.py')
    manifest = dict(kind='Geometry and mesh only', script_sha256=sha256(Path(__file__)),
                    checks=checks, source_cells=source['cells'], extended_cells=mesh['cells'],
                    source_nodes=len(source['xy']), extended_nodes=len(mesh['xy']),
                    annulus_method=provenance['mesher'], canvas_inches=[7.05,3.3], raster_dpi=600,
                    zoom_limits=zoom, visible_zoom_edges=dict(retained=len(retained), added=len(annulus),
                                                             shared=len(shared)),
                    software=dict(matplotlib=matplotlib.__version__, numpy=np.__version__),
                    inputs=[dict(path=str(path),sha256=sha256(path)) for path in files])
    (out/'source-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
