# Initial evidence audit

This audit reads existing saved outputs; it is not a new Fluent calculation.

- Cylinder: three mesh levels, 24,576 / 55,296 / 98,304 cells at Re=20.
  Raw force reports agree with the recorded drag coefficients and the
  0.01543057% medium-to-fine change (fine value used as denominator).
- Coarse-cylinder packaged transcript contains mesh checking only. The original
  continuation log records convergence at iteration 56 with all three monitored
  residuals below 1e-7. Its numerical excerpt is now included alongside that
  transcript; the complete source log stays outside Git.
- Medium/fine cylinder transcripts record convergence at iteration 51; both
  airfoil transcripts record convergence, at iterations 79 and 106.
- Final residual values and relative mass imbalances have been recalculated
  from the logs/reports in baseline-audit.json. All five final residual triplets
  meet 1e-7. Relative mass imbalance is small, but is not an accuracy bound.
- Fine-cylinder minimum orthogonal quality is 0.999985. Both airfoil meshes have
  minimum quality 0.0658879 and maximum aspect ratio about 80.84 at the trailing
  edge. The 5-degree airfoil transcript also records two early outlet reversed-flow
  messages; later convergence alone does not settle boundary sensitivity.
- Original scripts report final forces but do not configure force histories.
  Iterative force stability beyond residual convergence remains untested.
- The original portfolio's 44-entry content manifest was checked against the
  preserved import. New repository content has a separately regenerated manifest.

## Open questions

No formal grid error estimate, domain-size sensitivity, transient stability test,
matched benchmark comparison or aircraft validation is established. The airfoil
needs mesh repair and force sensitivity checks. Clean reproduction of the
portable script package remains pending. Concorde/ornithopter history was read
from the portfolio records; their live project states were not re-audited here.

## Operational readiness

All requested source folders and the D: Fluent executable exist. Imports and
Python are being prepared on E: with originals preserved. Consult the private
copy receipt and setup record for completed verification counts. The licence
utility returned `ERROR: No Feature List Retrieved` despite exit code zero;
Student licence checkout remains unconfirmed. No Fluent session was started.

Initial snapshot: about 490 GiB free on E:, 6.9 GiB available RAM, 8 physical / 16
logical CPU cores, no active Fluent/cortex process. These values are time-specific;
repeat checks immediately before substantial runs. Start with two solver
processes and only one session.

GitHub connector search and accessible repository listing found no
`Ocin12345/cfd-learning-lab`. Local preparation is not publication.


## Subsequent status

The setup copy receipt verified 12,715 files (443,409,577 bytes of file contents),
including the Python runtime and dependencies, with source/destination SHA-256
agreement and sources unchanged. Sources remain in place. Python 3.12.14 and
PyFluent 0.43.0 import from E:, and the project temp path resolves on E:.
Live Fluent Student launches subsequently succeeded for Lesson 1 and the first
domain-size case. The initial no-feature-list utility response remains preserved
as a historical diagnostic; it did not prevent those actual solver sessions.
See the lesson records for new evidence beyond the initial audit.
