#!/usr/bin/env python3
"""Extension scaffold: heteroskedastic outcomes.

This script operationalizes the professor's suggestion / Khan & Ugander perspective:
    k(e) = sigma1^2(e)/e + sigma0^2(e)/(1-e)
and compares propensity-only trimming with variance-aware trimming.

The simulation deliberately makes conditional variances functions of the propensity score e.
That is a transparent one-dimensional test bed; a fuller extension can let variances depend
on additional X not summarized by e(X).
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(__file__).resolve().parent / "outputs"
OUT.mkdir(exist_ok=True)

BETA_PAIRS = [(0.5,0.5),(0.5,1),(0.5,2),(0.5,4),(1,1),(1,2),(1,4),(2,2),(2,4),(4,4)]


def variance_scenarios(e):
    """Return named (sigma1^2, sigma0^2) scenarios."""
    one = np.ones_like(e)
    return {
        "homoskedastic": (one, one),
        # Same variance in both arms but much noisier at poor-overlap extremes.
        "symmetric_extreme_noise": (1 + 12*(e-0.5)**2, 1 + 12*(e-0.5)**2),
        # Treatment arm gets noisier as treatment becomes more common; control stays stable.
        # This generally makes the optimal k-threshold set asymmetric in propensity.
        "treatment_asymmetric_noise": (1 + 6*e**2, one),
        # A high-variance region near the center: shows that variance-aware trimming need not
        # simply discard the most extreme propensities.
        "central_noise_bump": (1 + 8*np.exp(-((e-0.5)/0.10)**2),
                               1 + 8*np.exp(-((e-0.5)/0.10)**2)),
    }


def objective(k, keep):
    p = keep.mean()
    if p == 0:
        return np.inf
    return np.mean(k * keep) / p**2


def best_symmetric_trim(e, k):
    # e is sorted because it was generated from increasing Beta quantiles.
    # Use prefix sums instead of recomputing a boolean mask for every alpha.
    alphas = np.linspace(0.0, 0.30, 3001)
    prefix = np.concatenate([[0.0], np.cumsum(k)])
    n = len(e)
    left = np.searchsorted(e, alphas, side="left")
    right = np.searchsorted(e, 1-alphas, side="right")
    count = right-left
    sumk = prefix[right]-prefix[left]
    vals = np.full_like(alphas, np.inf, dtype=float)
    ok = count > 0
    vals[ok] = n*sumk[ok]/(count[ok]**2)
    j = int(np.argmin(vals))
    a = float(alphas[j])
    keep = (e >= a) & (e <= 1-a)
    return a, float(vals[j]), keep


def best_k_threshold(k):
    """Minimize empirical E[k 1{k<=gamma}]/P(k<=gamma)^2 efficiently."""
    order = np.argsort(k)
    ks = k[order]
    cum = np.cumsum(ks)
    n = len(k)
    m = np.arange(1, n+1)
    vals = (cum/n) / (m/n)**2
    # Avoid pathological tiny target populations: theoretical optimum is interior here,
    # but requiring >=1% also makes finite-grid output stable.
    valid = m >= max(50, int(0.01*n))
    jj = np.argmin(np.where(valid, vals, np.inf))
    gamma = ks[jj]
    keep = k <= gamma
    return float(gamma), float(objective(k, keep)), keep


def summarize_rule(label, e, k, keep, v_full, v_fixed, v_kopt):
    kept_e = e[keep]
    v = objective(k, keep)
    return {
        "rule": label,
        "variance": v,
        "variance_ratio_to_full": v/v_full,
        "variance_ratio_to_fixed_0.1": v/v_fixed,
        "variance_ratio_to_variance_aware_opt": v/v_kopt,
        "retained_fraction": keep.mean(),
        "retained_e_min": kept_e.min() if keep.any() else np.nan,
        "retained_e_max": kept_e.max() if keep.any() else np.nan,
    }


def run_pair(beta, gamma, n_grid=200_000):
    # Deterministic midpoint quantiles: integration under Beta distribution without MC seed noise.
    q = (np.arange(n_grid) + 0.5) / n_grid
    e = stats.beta.ppf(q, beta, gamma)
    # Avoid literal 0/1 from any numerical underflow.
    e = np.clip(e, 1e-12, 1-1e-12)

    rows=[]
    for scen,(s1,s0) in variance_scenarios(e).items():
        k = s1/e + s0/(1-e)
        full = np.ones_like(e, dtype=bool)
        v_full = objective(k, full)

        keep_01 = (e>=0.1)&(e<=0.9)
        v_fixed = objective(k, keep_01)
        a_star, v_sym, keep_sym = best_symmetric_trim(e,k)
        g_star, v_k, keep_k = best_k_threshold(k)

        rules = [
            ("full", full),
            ("fixed_ps_0.1", keep_01),
            (f"optimal_symmetric_ps_alpha={a_star:.4f}", keep_sym),
            (f"optimal_variance_aware_k_gamma={g_star:.4g}", keep_k),
        ]
        for label, keep in rules:
            r=summarize_rule(label,e,k,keep,v_full,v_fixed,v_k)
            r.update({"beta":beta,"gamma":gamma,"scenario":scen,
                      "alpha_star_symmetric":a_star,"gamma_star_k":g_star})
            rows.append(r)
    return rows


def main():
    rows=[]
    for b,g in BETA_PAIRS:
        rows.extend(run_pair(b,g))
    df=pd.DataFrame(rows)
    df.to_csv(OUT/"extension_heteroskedasticity.csv",index=False)

    # A compact presentation focusing on the two rules you are most likely to discuss.
    compact=df[df["rule"].str.startswith(("fixed_ps","optimal_variance"))].copy()
    print(compact[["beta","gamma","scenario","rule","variance_ratio_to_fixed_0.1",
                   "variance_ratio_to_variance_aware_opt","retained_fraction",
                   "retained_e_min","retained_e_max"]].round(4).to_string(index=False))
    print("\nInterpretation: if the variance-aware retained e-range becomes asymmetric or excludes")
    print("an interior region, heteroskedasticity is doing something propensity-only trimming cannot capture.")

if __name__ == "__main__":
    main()
