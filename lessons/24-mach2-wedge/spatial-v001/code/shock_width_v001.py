"""Connected native-column10-90% transition width, with no extrapolation."""
import numpy as np
def connected_width(distance,normalized_pressure):
    d=np.asarray(distance,dtype=float);q=np.asarray(normalized_pressure,dtype=float)
    if d.ndim!=1 or q.shape!=d.shape or len(d)<2 or not np.isfinite(d).all() or not np.isfinite(q).all():raise ValueError("Finite equal-length profile arrays required")
    order=np.argsort(d);d,q=d[order],q[order]
    if np.any(np.diff(d)<=0):raise ValueError("Distinct profile coordinates required")
    candidates=np.flatnonzero((q[:-1]<=.5)&(q[1:]>.5))
    if not len(candidates):raise ValueError("No compression half-rise crossing")
    i=int(candidates[np.argmax((np.diff(q)/np.diff(d))[candidates])])
    left=i
    while left>=0 and q[left]>.1:
        if q[left]>.5:raise ValueError("10% bracket is separated by another half-rise")
        left-=1
    right=i+1
    while right<len(q) and q[right]<.9:
        if q[right]<.5:raise ValueError("90% bracket is separated by another half-rise")
        right+=1
    if left<0 or right>=len(q):raise ValueError("Missing connected10/90% bracket")
    def crossing(k,target):
        if not q[k]<=target<=q[k+1] or q[k+1]<=q[k]:raise ValueError("Invalid compression bracket")
        w=(target-q[k])/(q[k+1]-q[k]);return float(d[k]+w*(d[k+1]-d[k]))
    low,high,mid=crossing(left,.1),crossing(right-1,.9),crossing(i,.5)
    spacing=float(np.median(np.diff(d[max(0,left-1):min(len(d),right+2)])))
    if high<=low or spacing<=0:raise ValueError("Nonpositive transition width")
    return dict(lower_crossing_m=low,upper_crossing_m=high,midpoint_m=mid,width_m=high-low,projected_column_centroid_interval_m=spacing,width_in_projected_column_intervals=(high-low)/spacing,lower_bracket_indices=order[[left,left+1]].tolist(),upper_bracket_indices=order[[right-1,right]].tolist(),mid_bracket_indices=order[[i,i+1]].tolist())
