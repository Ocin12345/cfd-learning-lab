from pathlib import Path
import h5py,numpy as np,re
from plot_nozzle_native_v004 import indexed

def native_fluent_time(path):
    with h5py.File(path) as h:
        raw = h['settings/Data Variables'][()][0].decode()

    def val(key):
        a = re.findall('^\\(' + re.escape(key) + '\\s+([-+0-9.eE]+)\\)\\s*$', raw, re.M)
        assert len(a) == 1, (path, key, a)
        return float(a[0])
    return dict(time_s=val('flow-time'), step=val('time-step'), dt_s=val('physical-time-step'))

def fluent_wall_flux(case, data, report, names):
    raw = report.read_text()
    result = {}
    with h5py.File(case) as c, h5py.File(data) as d:
        mesh = c['meshes/1']
        z = mesh['faces/zoneTopology']
        nf = int(max(z['maxId']))
        xyz = mesh['nodes/coords']
        nn = max((int(np.asarray(a.attrs['maxId']).ravel()[0]) for a in xyz.values()))
        nodes = indexed(xyz, nn, 2)
        faces = np.full((nf, 2), -1, dtype=int)
        for group in mesh['faces/nodes'].values():
            a = group['nnodes']
            assert np.all(a[:] == 2)
            start = int(np.asarray(a.attrs.get('minId', group.attrs.get('minId', [1]))).ravel()[0]) - 1
            v = group['nodes'][:].reshape(-1, 2).astype(int) - 1
            assert np.all(faces[start:start + len(v)] == -1)
            faces[start:start + len(v)] = v
        assert np.all(faces >= 0) and np.all(faces < nn)
        lengths = np.linalg.norm(nodes[faces[:, 1]] - nodes[faces[:, 0]], axis=1)
        assert np.all(lengths > 0)
        q = indexed(d['results/1/phase-1/faces/SV_HEAT_FLUX'], nf, partial=True)
        mass = indexed(d['results/1/phase-1/faces/SV_FLUX'], nf)
        zn = z['name'][0].decode().split(';')
        assert set(names).issubset(zn)
        for name, start, end in zip(zn, z['minId'], z['maxId']):
            if name not in names:
                continue
            sl = slice(int(start) - 1, int(end))
            v = q[sl]
            m = mass[sl]
            assert np.isfinite(v).all()
            reported = re.findall('^\\s*' + re.escape(name) + '\\s+([-+0-9.eE]+)\\s*$', raw, re.M)
            assert len(reported) == 1
            Q = float(v.sum())
            assert np.isclose(Q, float(reported[0]), rtol=1e-07, atol=1e-10), (name, Q, reported)
            result[name] = dict(face_count=len(v), wall_length_m=float(lengths[sl].sum()), heat_per_depth_W_m=Q, report_heat_per_depth_W_m=float(reported[0]), absolute_mass_flux_per_depth_kg_ms=float(abs(m).sum()), net_mass_flux_per_depth_kg_ms=float(m.sum()), native_heat_representation='per-face heat rate, independently matched to Fluent Flux Report')
    return (result, [case, data, report])
