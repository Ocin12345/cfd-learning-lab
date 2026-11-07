"""Shared Mach 2 compression-corner geometry and weak-shock reference."""
from math import asin, atan, cos, degrees, pi, radians, sin, sqrt, tan

NX, NY = 200, 112
L, H, CORNER, DEPTH = 1., .82, .20, .01
THETA_DEG = 10.
THETA = radians(THETA_DEG)
GAMMA, R = 1.4, 287.05
CP = R * GAMMA / (GAMMA - 1.)
MW = 8314.46261815324 / R
MACH, P, T = 2., 101325., 300.
U = MACH * sqrt(GAMMA * R * T)
T0 = T * (1 + .5 * (GAMMA - 1.) * MACH**2)
P0 = P * (T0 / T)**(GAMMA / (GAMMA - 1.))
LIMITS = dict(mass_imbalance=.001, shock_angle_error_deg=1.5, pressure_ratio_error=.02,
              minimum_orthogonal_quality=.2, field_stability_relative=.002)
RESIDUALS = {"continuity": 1e-5, "x-velocity": 1e-5, "y-velocity": 1e-5, "energy": 1e-6}


def reference():
    """Solve the weak root of the theta–beta–Mach relation and normal jumps."""
    def turn(beta):
        m_normal_sq = (MACH * sin(beta))**2
        return atan(2 * (m_normal_sq - 1) /
                    (tan(beta) * (MACH**2 * (GAMMA + cos(2 * beta)) + 2)))
    a, b = asin(1. / MACH), pi / 2
    for k in range(80):
        low, high = (2 * a + b) / 3, (a + 2 * b) / 3
        if turn(low) < turn(high):
            a = low
        else:
            b = high
    a, b = asin(1. / MACH), (a + b) / 2
    if turn(b) < THETA:
        raise ValueError("The requested deflection has no attached shock")
    for k in range(70):
        mid = (a + b) / 2
        if turn(mid) < THETA:
            a = mid
        else:
            b = mid
    beta = (a + b) / 2
    mn2 = (MACH * sin(beta))**2
    p2p1 = 1 + 2 * GAMMA / (GAMMA + 1) * (mn2 - 1)
    rho2rho1 = (GAMMA + 1) * mn2 / ((GAMMA - 1) * mn2 + 2)
    m2n = sqrt((1 + .5 * (GAMMA - 1) * mn2) / (GAMMA * mn2 - .5 * (GAMMA - 1)))
    return dict(shock_angle_deg=degrees(beta), pressure_ratio=p2p1, density_ratio=rho2rho1,
                temperature_ratio=p2p1 / rho2rho1, downstream_mach=m2n / sin(beta - THETA))

