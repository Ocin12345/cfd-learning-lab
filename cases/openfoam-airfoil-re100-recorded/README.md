# Recorded airfoil dictionaries

These are the dictionaries actually used on the medium 383,562-cell airfoil.
Mesh, processor states, raw logs and saved fields remain outside Git. The
initial uniform-field run saved 2,000 iterations; U relaxation changed 0.7 to
0.9 for the continuation to 2,500, then remained 0.9 through 5,000. The control
dictionaries describe different stages, not interchangeable fresh-start recipes.
Preserve existing processor states rather than re-decomposing over them. See
[airfoil status](../../docs/AIRFOIL_MODEL_STATUS.md) for checks and limitations.
