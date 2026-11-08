import re,numpy as np
from plot_nozzle_native_v004 import foam_text,foam_list
from cavity_spatial_common_v001 import DT,DEPTH

def native_foam_time(path, t):
    raw = path.read_text()

    def val(key):
        a = re.findall('^\\s*' + key + '\\s+([-+0-9.eE]+);', raw, re.M)
        assert len(a) == 1, (path, key)
        return float(a[0])
    d = dict(time_s=val('value'), step=val('index'), dt_s=val('deltaT'))
    assert abs(d['time_s'] - t) < 1e-08 and abs(d['step'] * DT - t) < 1e-08 and (abs(d['dt_s'] - DT) < 1e-12), d
    return d

def scalar_patch(path, name, n, dimensions):
    raw = foam_text(path)
    d = re.search('dimensions\\s*\\[([^]]+)\\]', raw)
    assert d and np.array_equal(np.fromstring(d[1], sep=' '), dimensions), (path, 'dimensions')
    m = re.search('\\b' + re.escape(name) + '\\s*\\{([^{}]*)\\}', raw, re.S)
    assert m, (path, name)
    body = m[1]
    v = re.search('\\bvalue\\s+nonuniform\\s+List<scalar>\\s+(\\d+)\\s*\\((.*?)\\)\\s*;', body, re.S)
    if v:
        assert int(v[1]) == n
        values = np.fromstring(v[2], sep=' ')
    else:
        v = re.search('\\bvalue\\s+uniform\\s+([-+0-9.eE]+)\\s*;', body)
        assert v, (path, name)
        values = np.full(n, float(v[1]))
    assert values.shape == (n,) and np.isfinite(values).all()
    return values

def foam_walls(O, names, t):
    mesh = O / 'constant/polyMesh'
    if not mesh.exists():
        mesh = O / '0/polyMesh'
    n, body = foam_list(mesh / 'points')
    points = np.fromstring(body.replace('(', ' ').replace(')', ' '), sep=' ').reshape(n, 3)
    nf, body = foam_list(mesh / 'faces')
    matches = re.findall('(\\d+)\\s*\\(([^)]*)\\)', body)
    assert len(matches) == nf
    area = []
    for count, ids in matches:
        ids = np.fromstring(ids, sep=' ', dtype=int)
        assert len(ids) == int(count)
        assert np.all(ids >= 0) and np.all(ids < n)
        v = points[ids]
        area.append(float(np.linalg.norm(np.cross(v, np.roll(v, -1, axis=0)).sum(axis=0)) * 0.5))
    area = np.array(area)
    assert np.all(area > 0)
    boundary = foam_text(mesh / 'boundary')
    p = O / f'{t:g}'
    native_foam_time(p / 'uniform/time', t)
    report = O / 'wall-heat-postprocess.log'
    raw = report.read_text()
    assert re.search('^End\\s*$', raw, re.M)
    result = {}
    for name in names:
        block = re.search('\\b' + re.escape(name) + '\\s*\\{([^{}]*)\\}', boundary, re.S)
        assert block
        count = int(re.search('nFaces\\s+(\\d+)', block[1])[1])
        start = int(re.search('startFace\\s+(\\d+)', block[1])[1])
        a = area[start:start + count]
        assert len(a) == count
        q = scalar_patch(p / 'wallHeatFlux', name, count, [1, 0, -3, 0, 0, 0, 0])
        flux = scalar_patch(p / 'phi', name, count, [1, 0, -1, 0, 0, 0, 0])
        log = re.findall('for patch ' + re.escape(name) + ' = ([-+0-9.eE]+), ([-+0-9.eE]+), ([-+0-9.eE]+), ([-+0-9.eE]+)', raw)
        assert len(log) == 1
        Q = float(q @ a)
        reported = float(log[0][2])
        assert np.isclose(Q, reported, rtol=1e-09, atol=1e-13), (name, Q, reported)
        assert np.isclose(Q / a.sum(), float(log[0][3]), rtol=1e-09, atol=1e-10)
        result[name] = dict(face_count=count, wall_length_m=float(a.sum() / DEPTH), heat_per_depth_W_m=Q / DEPTH, report_heat_per_depth_W_m=reported / DEPTH, absolute_mass_flux_per_depth_kg_ms=float(abs(flux).sum() / DEPTH), net_mass_flux_per_depth_kg_ms=float(flux.sum() / DEPTH), native_heat_representation='W/m2 integrated over native polygon areas and divided by extrusion depth')
    return (result, [mesh / 'points', mesh / 'faces', mesh / 'boundary', p / 'wallHeatFlux', p / 'phi', p / 'uniform/time', report])
