"""Constant-property thermal settings for the existing planar slot jet."""
from jet_common_v001 import RHO, W

CP = 1006.
K = .0262
ALPHA = K / (RHO * CP)
COLD = 300.
HOT = 320.
DT = .005
LAST = 600
STRIDE = 20
ENERGY_RESIDUAL = 1e-8
ENERGY_BALANCE = .005

