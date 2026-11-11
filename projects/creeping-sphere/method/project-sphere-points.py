"""Project only the native sphere patch vertices to their prescribed radius."""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
POINT = re.compile(r"^\(("+NUMBER+r") ("+NUMBER+r") ("+NUMBER+r")\)$", re.M)
FACE = re.compile(r"^4\((\d+) (\d+) (\d+) (\d+)\)$", re.M)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def project(case, write):
    mesh = case / "constant/polyMesh"
    boundary = (mesh / "boundary").read_text(encoding="ascii")
    paths = {name: mesh / name for name in ("points", "faces", "owner", "neighbour", "boundary")}
    topology = {name: digest(path.read_bytes()) for name, path in paths.items() if name != "points"}
    patch = re.search(r"\bsphere\s*\{([^}]+)\}", boundary, re.S)[1]
    assert re.search(r"\btype\s+wall;", patch)
    count = int(re.search(r"nFaces\s+(\d+);", patch)[1])
    start = int(re.search(r"startFace\s+(\d+);", patch)[1])
    assert count == 1536
    faces = FACE.findall(paths["faces"].read_text(encoding="ascii"))
    assert len(faces) == 112128 and start+count == len(faces)
    vertices = {int(v) for face in faces[start:start+count] for v in face}
    assert len(vertices) == 1538
    source = paths["points"].read_bytes()
    points = source.decode("ascii")
    old = list(POINT.finditer(points))
    assert len(old) == 38450 and max(vertices) < len(old)

    chunks = []
    end = 0
    before = []
    after = []
    for index in sorted(vertices):
        point = old[index]
        coord = tuple(map(float, point.groups()))
        radius = math.hypot(*coord)
        assert abs(radius-1.) < .001, (index, radius)
        revised = tuple(format(value/radius, ".17g") for value in coord)
        chunks.extend((points[end:point.start()], "(" + " ".join(revised) + ")"))
        end = point.end()
        before.append(abs(radius-1.))
        after.append(abs(math.hypot(*map(float, revised))-1.))
    chunks.append(points[end:])
    revised_text = "".join(chunks)
    new = list(POINT.finditer(revised_text))
    assert len(new) == len(old)
    assert all(a[0] == b[0] for i, (a, b) in enumerate(zip(old, new)) if i not in vertices)
    assert max(after) <= 1e-14
    updated = revised_text.encode("ascii")
    receipt = {
        "sphere_vertices": len(vertices), "other_vertices_unchanged": len(old)-len(vertices),
        "radius_m": 1., "maximum_radius_error_before_m": max(before),
        "maximum_radius_error_after_m": max(after), "original_points_sha256": digest(source),
        "revised_points_sha256": digest(updated), "unchanged_topology_sha256": topology,
    }
    if write:
        with (case/"preprojection-points").open("xb") as stream:
            stream.write(source)
        temporary = paths["points"].with_name("points.projecting")
        with temporary.open("xb") as stream:
            stream.write(updated)
        assert paths["points"].read_bytes() == source
        temporary.replace(paths["points"])
        assert digest(paths["points"].read_bytes()) == receipt["revised_points_sha256"]
        assert all(digest(paths[name].read_bytes()) == value for name, value in topology.items())
        with (case/"sphere-projection.json").open("x", encoding="ascii") as stream:
            json.dump(receipt, stream, indent=2)
            stream.write("\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    project(args.case, args.write)
