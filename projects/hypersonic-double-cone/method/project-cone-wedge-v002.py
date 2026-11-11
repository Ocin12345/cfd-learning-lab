"""Restore the intended planar wedge coordinates without changing x, y or topology."""
import hashlib
import json
import math
import re
from pathlib import Path

mesh = Path("constant/polyMesh")
points = mesh/"points"
old = points.read_bytes()
pattern = re.compile(rb"^\(([-+0-9.eE]+) ([-+0-9.eE]+) ([-+0-9.eE]+)\)$", re.M)
changes = []
axis_changes = []
vertex_index = -1
faces = [tuple(map(int, face.split())) for face in re.findall(rb"^\d+\(([\d ]+)\)$",(mesh/"faces").read_bytes(),re.M)]
assert len(faces)==129795
left = {v for f in faces[65025:97410] for v in f}
right = {v for f in faces[97410:] for v in f}
axis = left & right
assert len(axis)==128 and len(left | right)==65408
half_angle = math.tan(math.radians(.1))


def fix_point(match):
    global vertex_index
    vertex_index += 1
    x, y, z = match.groups()
    radius, thickness = float(y), float(z)
    assert radius >= 0
    if vertex_index in axis:
        if radius or thickness:
            assert float(x) < 0
            axis_changes.append(math.hypot(radius,thickness))
        return b"("+x+b" 0 0)"
    sign = -1 if vertex_index in left else 1
    assert radius > 0 and sign*thickness > 0
    corrected = sign*radius*half_angle
    changes.append(abs(corrected-thickness))
    return b"("+x+b" "+y+b" "+format(corrected,".17g").encode()+b")"


new, n = pattern.subn(fix_point,old)
assert n == 65408 and 0 < max(changes) < 5e-8
assert len(axis_changes)==2 and max(axis_changes)<5e-6
with Path("preprojection-wedge-points").open("xb") as f:
    f.write(old)
temp = mesh/"points.projecting"
with temp.open("xb") as f:
    f.write(new)
assert points.read_bytes() == old
temp.replace(points)
assert points.read_bytes() == new
record = {"native_point_count":n, "max_azimuthal_coordinate_change_m":max(changes),
          "half_angle_degrees":.1, "off_axis_shared_nodes":len(axis_changes),
          "max_axis_node_change_m":max(axis_changes), "preprojection_sha256":hashlib.sha256(old).hexdigest(),
          "points_sha256":hashlib.sha256(new).hexdigest()}
print(json.dumps(record,indent=2))
