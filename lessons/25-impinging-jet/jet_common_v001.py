"""Parameters shared by both planar slot-jet runs."""
import numpy as np

W = .02
H = 3 * W
HALF = 6 * W
NX, NY, NSLOT = 240, 60, 20
DX = W / NSLOT
DEPTH = .01
RHO, MU = 1.225, 1.7894e-5
RE = 100.
UMEAN = RE * MU / (RHO * W)
Q_REF = .5 * RHO * UMEAN**2
X = (np.arange(NX + 1) - NX / 2) * DX
Y = np.arange(NY + 1) * DX
X_SLOT = (np.arange(NSLOT) - (NSLOT - 1) / 2) * DX
V_SLOT = -1.5 * UMEAN * (1 - (2 * X_SLOT / W)**2)
INLET_MASS = -RHO * DX * float(V_SLOT.sum())
RESIDUAL = 1e-6
MASS_BALANCE = 1e-4
FIELD_STABILITY = .002
FLOW_SPLIT_ERROR = .005
ZONES = ("inlet", "left", "right", "plate", "ceiling")

