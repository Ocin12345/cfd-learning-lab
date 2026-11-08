import math
def compare_restart_settings(actual, base):
    roundoff = []
    roundoff = []

    def same(a, b, path):
        if isinstance(a, dict) and isinstance(b, dict):
            assert a.keys() == b.keys(), (path, a.keys(), b.keys())
            for k in a:
                same(a[k], b[k], path + '/' + k)
        elif isinstance(a, list) and isinstance(b, list):
            assert len(a) == len(b), path
            for i, (x, y) in enumerate(zip(a, b)):
                same(x, y, path + '/' + str(i))
        elif isinstance(a, (int, float)) and (not isinstance(a, bool)) and isinstance(b, (int, float)) and (not isinstance(b, bool)):
            assert math.isclose(a, b, rel_tol=1e-14, abs_tol=1e-30), (path, a, b)
            if a != b:
                roundoff.append(dict(path=path, actual=a, baseline=b))
        else:
            assert a == b, (path, a, b)
    for key in ['general', 'energy', 'viscous', 'material', 'boundaries', 'references', 'methods', 'controls']:
        same(actual[key], base[key], key)
    assert actual['methods']['p_v_coupling']['flow_scheme'] == 'PISO'
    return dict(passed=True, numeric_relative_tolerance=1e-14, roundoff_only=roundoff, compared_keys=['general', 'energy', 'viscous', 'material', 'boundaries', 'references', 'methods', 'controls'])
