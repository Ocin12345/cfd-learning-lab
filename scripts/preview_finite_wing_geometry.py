"""Plot actual exported wing mesh/profile geometry only; no CFD fields implied."""
from pathlib import Path
import argparse
import importlib.util
import json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out, geometry_path = args.out.resolve(), args.geometry.resolve()
    if out.drive.upper() != 'E:' or out.exists():
        raise ValueError('Use a fresh E: geometry-preview directory')
    reader = Path(__file__).with_name('plot_finite_wing_native.py')
    spec = importlib.util.spec_from_file_location('wing_native', reader)
    native = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(native)
    np, plt = native.np, native.plt
    geometry = native.TetGeometry(geometry_path)
    with np.load(geometry_path, allow_pickle=False) as data:
        profile = data['profile_xyz_incidence'].copy()
        raw_profile = data['profile_xy']
    leading, trailing = profile[np.argmin(raw_profile[:, 0])], profile[np.argmax(raw_profile[:, 0])]
    chord = trailing - leading
    angle = float(np.degrees(np.arctan2(chord[1], chord[0])))
    chord_length = float(np.linalg.norm(chord))
    quarter = leading + .25 * chord
    span = [float(geometry.wing_vertices[:, :, 2].min()), float(geometry.wing_vertices[:, :, 2].max())]
    native.require(abs(angle + 8) < 1e-8 and abs(chord_length - 1) < 1e-8
                   and np.max(np.abs(quarter - [.25, 0, 0])) < 1e-8 and np.allclose(span, [-2, 2]),
                   'Exported profile/full-span geometry differs from the approved coordinate contract')
    out.mkdir(parents=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig = plt.figure(figsize=(12.8, 6.3))
    ax = fig.add_axes([.015, .13, .61, .74], projection='3d')
    ax.add_collection3d(native.Poly3DCollection(native.display(geometry.wing_vertices),
        facecolors='#c4d6e2', edgecolors=(.2, .35, .46, .20), linewidths=.12))
    for z in span:
        outline = profile.copy()
        outline[:, 2] = z
        points = native.display(np.vstack([outline, outline[0]]))
        ax.plot(points[:, 0], points[:, 1], points[:, 2], color='#193f57', lw=1.0)
    native.decorate_3d(ax, 'Actual closed wing surface mesh')
    ax.set(xlim=(-.1, 1.2), ylim=(-2.3, 2.3), zlim=(-.38, .38))
    ax.set_box_aspect((1.3, 4.6, .76))
    ax.view_init(elev=25, azim=-47)
    side = fig.add_axes([.68, .31, .3, .40])
    closed = np.vstack([profile, profile[0]])
    side.plot(closed[:, 0], closed[:, 1], color='#1d5a79', lw=1.8)
    side.plot([leading[0], trailing[0]], [leading[1], trailing[1]], color='#87929a', ls='--', lw=.9)
    side.scatter([.25], [0], color='#c57032', s=26, zorder=5)
    side.annotate('Quarter-chord pivot', (.25, 0), xytext=(.31, .14), arrowprops={'arrowstyle': '-', 'color': '#9a734b'}, fontsize=9)
    side.set(xlim=(-.05, 1.05), ylim=(-.20, .22), aspect='equal', xlabel='$x$ (m)', ylabel='$y$ (m)')
    side.set_title('Exported profile at 8° clockwise pitch', fontsize=10, pad=14)
    side.spines[['top', 'right']].set_visible(False)
    side.grid(alpha=.15)
    fig.suptitle('Finite-span wing · meshed geometry only', fontsize=17, y=.975)
    fig.text(.5, .055, f'{len(geometry.tets):,} fluid tetrahedra · {len(geometry.points):,} nodes · {len(geometry.wing_faces):,} wing triangles',
              ha='center', fontsize=10, color='#4b5963')
    fig.text(.5, .018, 'Chord 1 m · full span 4 m · no computed pressure, velocity or vortex field shown', ha='center', fontsize=9, color='#64717a')
    fig.savefig(out / 'finite-wing-actual-mesh-and-profile.png', dpi=180)
    plt.close(fig)
    checks = dict(actual_exported_mesh=True, CFD_fields_shown=False, fluid_tetrahedra=len(geometry.tets),
        nodes=len(geometry.points), wing_triangles=len(geometry.wing_faces), chord_m=chord_length,
        measured_pitch_degrees=angle, quarter_chord_pivot_m=quarter.tolist(), full_span_bounds_m=span,
        minimum_tet_volume_m3=float(geometry.volumes.min()), source_geometry_sha256=native.sha(geometry_path))
    (out / 'checks.json').write_text(json.dumps(checks, indent=2) + '\n', encoding='utf-8')
    (out / 'executed-preview.py').write_bytes(Path(__file__).read_bytes())
    (out / 'executed-native-geometry-reader.py').write_bytes(reader.read_bytes())
    (out / 'input-hash-manifest.json').write_text(json.dumps(dict(inputs=[dict(path=str(geometry_path), sha256=native.sha(geometry_path))],
        preview_sha256=native.sha(Path(__file__)), reader_sha256=native.sha(reader)), indent=2) + '\n', encoding='utf-8')
    (out / 'CAPTION.md').write_text('Actual exported finite-wing surface mesh and incidence-transformed profile. '
        'The closed symmetric NACA0012-like profile has a rounded trailing cap; it is not asserted to be exact NACA0012. '
        'Clockwise 8 degree rotation is about the span axis through (0.25,0,0) m; chord is 1 m and span is -2..2 m. '
        'The surface mesh is drawn from actual NPZ connectivity. This is meshed geometry only, with no computed '
        'pressure, velocity, downwash or vortices represented. The fluid tetrahedra and wing-face ownership passed '
        'connectivity/positive-volume checks; solver-native quality checks remain separate.\n', encoding='utf-8')
    files = [p for p in out.iterdir() if p.is_file()]
    (out / 'output-hash-manifest.json').write_text(json.dumps([dict(path=p.name, bytes=p.stat().st_size, sha256=native.sha(p)) for p in files], indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(image=str(out / 'finite-wing-actual-mesh-and-profile.png'), **checks), indent=2))


if __name__ == '__main__':
    main()
