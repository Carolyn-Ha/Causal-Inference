#!/usr/bin/env python3
"""Replicate Crump et al. (2009), Table 1.

The paper studies E ~ Beta(beta, gamma) and, under homoskedasticity,
  V(alpha) / sigma^2 = E[ 1{alpha<=E<=1-alpha} * (1/E + 1/(1-E)) ]
                       / P(alpha<=E<=1-alpha)^2.
The factor sigma^2 cancels in all variance ratios.

Two outputs are intentionally reported:
  1) literal_formula: mathematically exact where finite.  The full-sample variance
     is infinite if beta<=1 or gamma<=1.
  2) paper_matching_diagnostic: a small endpoint floor (default 1e-4) is used
     only for the *full-sample* integral.  This reproduces the published finite
     numbers fairly closely in the divergent cases, but the paper does not state
     such a floor; therefore this is a diagnostic, not a claim about the authors'
     implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy import integrate, optimize, stats

OUT = Path(__file__).resolve().parent / "outputs"
OUT.mkdir(exist_ok=True)

# Published Table 1 cells (rows beta, columns gamma) for comparison.
PUBLISHED = {
    (0.5, 0.5): (13.38, 1.70, 1.00),
    (0.5, 1.0): (11.68, 1.64, 1.00),
    (0.5, 2.0): (13.71, 1.71, 1.00),
    (0.5, 4.0): (12.83, 1.58, 1.04),
    (1.0, 1.0): (2.68, 1.39, 1.00),
    (1.0, 2.0): (2.65, 1.39, 1.00),
    (1.0, 4.0): (3.36, 1.47, 1.01),
    (2.0, 2.0): (1.11, 1.09, 1.00),
    (2.0, 4.0): (1.16, 1.12, 1.00),
    (4.0, 4.0): (1.02, 1.02, 1.00),
}


def beta_mass(alpha: float, beta: float, gamma: float) -> float:
    """P(alpha <= E <= 1-alpha), E~Beta(beta,gamma)."""
    if not (0 <= alpha < 0.5):
        return 0.0
    return float(stats.beta.cdf(1 - alpha, beta, gamma) -
                 stats.beta.cdf(alpha, beta, gamma))


def variance_alpha(alpha: float, beta: float, gamma: float) -> float:
    """V_S,alpha / sigma^2 from Eq. in Sec. 3.3; finite for alpha>0."""
    if not (0 < alpha < 0.5):
        raise ValueError("variance_alpha expects 0 < alpha < 0.5")
    p = beta_mass(alpha, beta, gamma)
    if p <= 0:
        return math.inf

    def integrand(e: float) -> float:
        return (1.0 / e + 1.0 / (1.0 - e)) * stats.beta.pdf(e, beta, gamma)

    numerator, _ = integrate.quad(
        integrand, alpha, 1 - alpha,
        epsabs=1e-10, epsrel=1e-10, limit=500
    )
    return float(numerator / (p * p))


def full_variance_exact(beta: float, gamma: float) -> float:
    """V_S/sigma^2 under literal full support.

    For E~Beta(beta,gamma), E(1/E) is finite iff beta>1 and
    E(1/(1-E)) is finite iff gamma>1.
    """
    if beta <= 1 or gamma <= 1:
        return math.inf
    # E[1/E] = (beta+gamma-1)/(beta-1)
    # E[1/(1-E)] = (beta+gamma-1)/(gamma-1)
    common = beta + gamma - 1
    return common / (beta - 1) + common / (gamma - 1)


def full_variance_endpoint_floor(beta: float, gamma: float,
                                 eps: float = 1e-4,
                                 renormalize: bool = True) -> float:
    """Numerical endpoint-regularized version of the full-sample variance.

    This is NOT part of the published method.  It is provided because the
    published Table 1 reports finite full-sample values even for beta<=1,
    where the literal expectation diverges.

    renormalize=True interprets [eps,1-eps] as an implicit retained sample.
    """
    if not (0 < eps < 0.5):
        raise ValueError("eps must be in (0, 0.5)")

    def integrand(e: float) -> float:
        return (1.0 / e + 1.0 / (1.0 - e)) * stats.beta.pdf(e, beta, gamma)

    numerator, _ = integrate.quad(
        integrand, eps, 1 - eps,
        epsabs=1e-10, epsrel=1e-10, limit=500
    )
    if not renormalize:
        return float(numerator)
    p = float(stats.beta.cdf(1 - eps, beta, gamma) - stats.beta.cdf(eps, beta, gamma))
    return float(numerator / (p * p))


def optimal_alpha(beta: float, gamma: float) -> tuple[float, float]:
    """Numerically minimize V(alpha) over symmetric trimming rules."""
    # Small positive lower bound avoids divergent endpoints during optimization.
    res = optimize.minimize_scalar(
        lambda a: variance_alpha(a, beta, gamma),
        bounds=(1e-7, 0.499),
        method="bounded",
        options={"xatol": 1e-11, "maxiter": 1000},
    )
    if not res.success:
        raise RuntimeError(res.message)
    return float(res.x), float(res.fun)


def one_pair(beta: float, gamma: float, endpoint_eps: float = 1e-4) -> dict:
    a_star, v_star = optimal_alpha(beta, gamma)
    v01 = variance_alpha(0.01, beta, gamma)
    v10 = variance_alpha(0.10, beta, gamma)
    vfull = full_variance_exact(beta, gamma)
    vfull_diag = full_variance_endpoint_floor(beta, gamma, endpoint_eps, renormalize=True)
    pub = PUBLISHED[(beta, gamma)]
    return {
        "beta": beta,
        "gamma": gamma,
        "alpha_star": a_star,
        "Vfull_over_Vopt_literal": vfull / v_star if np.isfinite(vfull) else np.inf,
        "Vfull_over_Vopt_endpoint_floor": vfull_diag / v_star,
        "V001_over_Vopt": v01 / v_star,
        "V010_over_Vopt": v10 / v_star,
        "paper_Vfull_over_Vopt": pub[0],
        "paper_V001_over_Vopt": pub[1],
        "paper_V010_over_Vopt": pub[2],
        "absdiff_floor_full": abs(vfull_diag / v_star - pub[0]),
        "absdiff_001": abs(v01 / v_star - pub[1]),
        "absdiff_010": abs(v10 / v_star - pub[2]),
    }


def make_long_table(endpoint_eps: float = 1e-4) -> pd.DataFrame:
    pairs = list(PUBLISHED.keys())
    return pd.DataFrame([one_pair(b, g, endpoint_eps) for b, g in pairs])


def make_paper_layout(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    """Triangular beta x gamma layout matching the published table."""
    vals = [0.5, 1.0, 2.0, 4.0]
    out = pd.DataFrame(index=vals, columns=vals, dtype=float)
    for _, r in df.iterrows():
        out.loc[r["beta"], r["gamma"]] = r[value_col]
    out.index.name = "beta"
    out.columns.name = "gamma"
    return out


def main() -> None:
    endpoint_eps = 1e-4
    df = make_long_table(endpoint_eps)
    df.to_csv(OUT / "table1_replication_long.csv", index=False)

    # Save easy-to-read versions of each of the three published ratio panels.
    panels = {
        "full_literal": "Vfull_over_Vopt_literal",
        "full_endpoint_floor": "Vfull_over_Vopt_endpoint_floor",
        "trim_001": "V001_over_Vopt",
        "trim_010": "V010_over_Vopt",
    }
    for name, col in panels.items():
        make_paper_layout(df, col).round(4).to_csv(OUT / f"table1_{name}.csv")

    print("\nCrump et al. (2009) Table 1 replication\n")
    show_cols = [
        "beta", "gamma", "alpha_star",
        "Vfull_over_Vopt_literal", "Vfull_over_Vopt_endpoint_floor",
        "V001_over_Vopt", "V010_over_Vopt",
        "paper_Vfull_over_Vopt", "paper_V001_over_Vopt", "paper_V010_over_Vopt",
    ]
    with pd.option_context("display.max_columns", None, "display.width", 180):
        print(df[show_cols].round(4).to_string(index=False))

    print("\nImportant diagnostic:")
    print("  Under the literal Beta model, full-sample V is infinite whenever beta<=1 or gamma<=1.")
    print("  The published table nevertheless reports finite values in those cells.")
    print(f"  The endpoint-floor column uses eps={endpoint_eps:g} only as a reproduction diagnostic.")
    print("  Do not describe that floor as an author-specified method unless you find original code documenting it.")


if __name__ == "__main__":
    main()
