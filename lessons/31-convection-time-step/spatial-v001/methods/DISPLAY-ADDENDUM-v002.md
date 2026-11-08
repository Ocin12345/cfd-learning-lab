# Display-only correction, version 002

Both full solver continuations are accepted and closed. No CFD rerun is authorized by this display correction.

Version 001's streamline integrator seeded the near-stagnation centre and produced overlapping tiny arrowheads (a white dot). Version 002 masks streamline integration where the local native speed is below 0.5% of the shared fixed 10 mm/s colour limit and removes corresponding seeds. It leaves all native cell values and colour maps unchanged. Curves remain instantaneous streamlines, not particle paths. Endpoint native face-heat markers are unconnected; no two-point transient history is implied.

The original v001 render, metadata-only NumPy-scalar repair, images and complete MP4/GIF remain preserved. Version 002 uses strict finite JSON with an explicit NumPy scalar converter. All 164 collected native states, CSV calculations and source hashes are unchanged and will be compared before selection. Two display tests pass: centre masking without input mutation and quiet-field seed suppression. No solver acceptance threshold changes.
