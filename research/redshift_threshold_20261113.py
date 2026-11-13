"""Conditional truth-labelled redshift diagnostics, not a Euclid Q1 result."""
import math
def residuals(z, truth):
    if not z or len(z)!=len(truth) or any(not math.isfinite(v) for v in list(z)+list(truth)) or any(v<=-1 for v in truth):
        raise ValueError("invalid redshifts")
    return [abs(a-b)/(1+b) for a,b in zip(z,truth)]
def wilson(k,n,z=1.959963984540054):
    if n<1 or k<0 or k>n or z<=0: raise ValueError("invalid counts")
    p=k/n;d=1+z*z/n;c=(p+z*z/(2*n))/d
    h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return max(0.,c-h),min(1.,c+h)
def threshold_sweep(r,thresholds):
    if not r or not thresholds or any(not math.isfinite(x) or x<0 for x in r) or any(not math.isfinite(t) or t<=0 for t in thresholds):
        raise ValueError("invalid inputs")
    return [(t,sum(x>t for x in r),wilson(sum(x>t for x in r),len(r))) for t in thresholds]
def leave_one_source_out(r,ids,t):
    if len(r)!=len(ids) or len(r)<2 or len(set(ids))!=len(ids) or t<=0:
        raise ValueError("require unique sources")
    return {sid:sum(x>t for j,x in enumerate(r) if i!=j)/(len(r)-1) for i,sid in enumerate(ids)}
