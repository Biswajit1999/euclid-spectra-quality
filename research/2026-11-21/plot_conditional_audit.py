"""Generate a labelled synthetic uncertainty coverage figure; not Euclid data."""
from pathlib import Path
import importlib.util
import numpy as np
import matplotlib.pyplot as plt
spec = importlib.util.spec_from_file_location("audit", Path(__file__).with_name("conditional_audit.py"))
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
z_ref = np.linspace(.2,2,100)
z_est = z_ref + .06*np.sin(np.arange(100)*.8)
sigma = np.full(100,.05)
mult = np.linspace(.25,3,40)
coverage = [audit.interval_coverage(z_est,z_ref,sigma,x) for x in mult]
plt.plot(mult,coverage,label="Empirical synthetic coverage")
plt.xlabel("Quoted uncertainty multiplier")
plt.ylabel("Fraction covered")
plt.title("Synthetic redshift interval audit — not Euclid data")
plt.ylim(0,1.02)
plt.legend()
plt.tight_layout()
plt.savefig("synthetic_euclid_coverage_audit.png",dpi=160)
