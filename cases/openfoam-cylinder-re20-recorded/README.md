# Recorded OpenFOAM dictionaries

These are the small input dictionaries from the actual paired-cylinder run.
The mesh and native results remain outside Git. The first control dictionary
starts at zero and caps at 1,200 iterations; the second resumes from the saved
latest state to 3,000. Preserve that distinction when reproducing the continuation.
Use scripts/prepare_openfoam_cylinder.py to prepare a fresh external case,
then convert and check the mesh before solving. Do not interpret the steady
iteration coordinate as a physical flow-time history.
