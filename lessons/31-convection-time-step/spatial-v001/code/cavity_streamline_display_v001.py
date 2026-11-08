from pathlib import Path
import numpy as np

def streamline_display(U,centres,seed,threshold_m_s):
    """Mask tiny native speeds for drawing only; leave solver fields unchanged."""
    speed=np.linalg.norm(U,axis=-1)
    mask=speed<threshold_m_s
    ix=np.abs(centres[:,None]-seed[:,0]).argmin(axis=0)
    iy=np.abs(centres[:,None]-seed[:,1]).argmin(axis=0)
    return np.ma.array(U[...,0],mask=mask,copy=True),np.ma.array(U[...,1],mask=mask,copy=True),seed[~mask[iy,ix]].copy()
