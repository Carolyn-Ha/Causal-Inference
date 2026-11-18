"""Paired oracle-vs-measured network error extension to Cai stochastic seeding.

The true graph always generates assignment mechanisms and outcomes. Only the
analyst's one-hop probabilities are computed on edge-deleted graphs.

IMPORTANT: estimated oracle effects are NOT the true effect. An independent,
configurable Monte Carlo reference estimates the true stochastic-seeding effect.
Paired measured-minus-oracle differences quantify the incremental impact of
network measurement error without needing a reference truth estimate.

Run (from project root):
 python official_python_extension.py --replicates 100 --truth-replicates 1000 --output outputs/pilot
 python official_python_extension.py --replicates 5000 --truth-replicates 20000 --output outputs/full --resume
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from numba import njit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from official_python_replication import (
    load_author_villages, choose_seeds, adoption, ALPHAS, BETAS, DEGREE_EFFECTS,
)

RATES = (0.0, 0.1, 0.2, 0.3, 0.5)
METHODS = ('DM', 'HT', 'Hajek')
WEIGHT_TYPES = ('oracle', 'measured')
Z90 = 1.6448536269514722


@njit
def thinned_nomination_probability(A, n, rate, sa, sb):
    """Draw ONE corrupted directed multigraph and compute pair probability.

    Each directed edge occurrence is deleted independently with probability
    ``rate``. Individuals with zero remaining out-degree cannot nominate.
    A selected pair can have probability zero under the measured graph.
    """
    q = np.zeros(n, np.float64)
    n_valid = 0
    for i in range(n):
        d = 0
        row = np.zeros(n, np.int64)
        for j in range(n):
            for _ in range(A[i, j]):
                if np.random.random() >= rate:
                    row[j] += 1
                    d += 1
        if d > 0:
            n_valid += 1
            for j in range(n):
                if row[j] > 0:
                    q[j] += row[j] / d
    if n_valid == 0:
        return 0.0
    q /= n_valid
    norm = 1.0 - (q * q).sum()
    if norm <= 0.0:
        return 0.0
    return 2.0 * q[sa] * q[sb] / norm


@njit
def estimate_once(ys, ws, pa, pb):
    """DM, HT and Hajek under the mixed two-arm proposal 0.5 pa + 0.5 pb.

    Returns estimates and nominal standard errors. Nominal standard errors
    treat the supplied probability model as correct; they are NOT corrected
    for graph measurement error or model uncertainty.
    """
    N = len(ys)
    n1 = 0
    n0 = 0
    sum1 = 0.0
    sum0 = 0.0
    wa = np.empty(N)
    wb = np.empty(N)
    ht_scores = np.empty(N)
    dena = 0.0
    denb = 0.0
    numa = 0.0
    numb = 0.0
    for j in range(N):
        if ws[j]:
            n1 += 1
            sum1 += ys[j]
        else:
            n0 += 1
            sum0 += ys[j]
        d = 0.5 * (pa[j] + pb[j])
        wa[j] = pa[j] / d
        wb[j] = pb[j] / d
        ht_scores[j] = (wa[j] - wb[j]) * ys[j]
        dena += wa[j]
        denb += wb[j]
        numa += wa[j] * ys[j]
        numb += wb[j] * ys[j]
    mu1 = sum1 / n1 if n1 else np.nan
    mu0 = sum0 / n0 if n0 else np.nan
    dm = mu1 - mu0
    dm_var = np.nan
    if n1 > 1 and n0 > 1:
        v1 = 0.0
        v0 = 0.0
        for j in range(N):
            if ws[j]:
                v1 += (ys[j] - mu1) ** 2
            else:
                v0 += (ys[j] - mu0) ** 2
        dm_var = v1 / (n1 * (n1 - 1)) + v0 / (n0 * (n0 - 1))
    ht = ht_scores.mean()
    ht_var = 0.0
    for j in range(N):
        ht_var += (ht_scores[j] - ht) ** 2
    ht_var /= N * N
    mua = numa / dena if dena > 0 else np.nan
    mub = numb / denb if denb > 0 else np.nan
    haj = mua - mub
    haj_var = 0.0
    if np.isfinite(haj):
        for j in range(N):
            influence = wa[j] * (ys[j] - mua) - wb[j] * (ys[j] - mub)
            haj_var += influence * influence
        haj_var /= N * N
    else:
        haj_var = np.nan
    return np.array((dm, ht, haj)), np.sqrt(np.array((dm_var, ht_var, haj_var)))


@njit
def independent_reference(A, Din, Dout, Q, Ns, alpha, gamma, beta,
                          N, reps, seed):
    """Independent MC reference for E[Y(onehop)-Y(random)] with N villages.

    Each replicate randomly samples N of the 150 villages and independently
    draws two potential adoption outcomes for each selected village. This is
    the same superpopulation-over-villages simulation target as the original
    Python port's averaged per-replicate 'truth', but uses disjoint RNG runs.
    """
    np.random.seed(seed)
    M = len(Ns)
    differences = np.empty(reps)
    means_a = np.empty(reps)
    means_b = np.empty(reps)
    for r in range(reps):
        order = np.random.permutation(M)
        total_a = 0.0
        total_b = 0.0
        for j in range(N):
            vi = order[j]
            n = Ns[vi]
            sa, sb = choose_seeds(Q[vi], n, True)
            ra, rb = choose_seeds(Q[vi], n, False)
            total_a += adoption(A[vi], Dout[vi], Din[vi], n, sa, sb, alpha, gamma, beta)
            total_b += adoption(A[vi], Dout[vi], Din[vi], n, ra, rb, alpha, gamma, beta)
        means_a[r] = total_a / N
        means_b[r] = total_b / N
        differences[r] = (total_a - total_b) / N
    return differences, means_a, means_b


@njit
def run_cell(A, Din, Dout, Q, Norm, Ns, alpha, gamma, beta,
             N, reps, rates, seed):
    """Keep the whole observation fixed for paired oracle/measured comparison.

    estimates: [rep, rate, {oracle, measured}, {DM,HT,Hajek}, {value,SE}].
    Zero-probability counts are number of observed pairs assigned p_A=0 by
    the corrupted network, even if the oracle p_A is positive.
    """
    np.random.seed(seed)
    M = len(Ns)
    R = len(rates)
    estimates = np.empty((reps, R, 2, 3, 2), np.float64)
    zero_count = np.zeros((reps, R), np.int64)
    zero_count_oracle_positive = np.zeros((reps, R), np.int64)
    for r in range(reps):
        order = np.random.permutation(M)
        y = np.empty(N)
        z = np.empty(N, np.uint8)
        seed_a = np.empty(N, np.int64)
        seed_b = np.empty(N, np.int64)
        indices = np.empty(N, np.int64)
        pa_true = np.empty(N)
        pb = np.empty(N)
        for j in range(N):
            vi = order[j]
            indices[j] = vi
            n = Ns[vi]
            z[j] = np.random.randint(2)
            sa, sb = choose_seeds(Q[vi], n, True)
            ra, rb = choose_seeds(Q[vi], n, False)
            # Draw both potential outcomes as in the original replication.
            yn = adoption(A[vi], Dout[vi], Din[vi], n, sa, sb, alpha, gamma, beta)
            yr = adoption(A[vi], Dout[vi], Din[vi], n, ra, rb, alpha, gamma, beta)
            if z[j]:
                seed_a[j], seed_b[j] = sa, sb
                y[j] = yn
            else:
                seed_a[j], seed_b[j] = ra, rb
                y[j] = yr
            pa_true[j] = 2.0 * Q[vi, seed_a[j]] * Q[vi, seed_b[j]] / Norm[vi]
            pb[j] = 2.0 / (n * (n - 1))
        oracle_vals, oracle_ses = estimate_once(y, z, pa_true, pb)
        for k in range(R):
            rate = rates[k]
            pm = pa_true.copy()
            if rate > 0:
                for j in range(N):
                    vi = indices[j]
                    pm[j] = thinned_nomination_probability(
                        A[vi], Ns[vi], rate, seed_a[j], seed_b[j]
                    )
            for j in range(N):
                if pm[j] == 0:
                    zero_count[r, k] += 1
                    if pa_true[j] > 0:
                        zero_count_oracle_positive[r, k] += 1
            measured_vals, measured_ses = estimate_once(y, z, pm, pb)
            for m in range(3):
                estimates[r, k, 0, m, 0] = oracle_vals[m]
                estimates[r, k, 0, m, 1] = oracle_ses[m]
                estimates[r, k, 1, m, 0] = measured_vals[m]
                estimates[r, k, 1, m, 1] = measured_ses[m]
    return estimates, zero_count, zero_count_oracle_positive


def _stats(vals, ses, tau):
    valid = np.isfinite(vals)
    x = vals[valid]
    s = ses[valid]
    if not len(x):
        raise ValueError('All estimator values are non-finite')
    errors = x - tau
    ci_valid = np.isfinite(s)
    return {
        'mean_estimate': float(np.mean(x)),
        'bias': float(np.mean(errors)),
        'variance': float(np.var(x, ddof=1)) if len(x) > 1 else np.nan,
        'rmse': float(np.sqrt(np.mean(errors**2))),
        'mean_nominal_se': float(np.mean(s[ci_valid])) if np.any(ci_valid) else np.nan,
        'coverage_90': float(np.mean(np.abs(errors[ci_valid]) <= Z90*s[ci_valid])) if np.any(ci_valid) else np.nan,
        'n_valid': int(len(x)),
        'n_ci_valid': int(np.sum(ci_valid)),
    }


def summarize_cell(estimates, zero, zero_positive, tau, tau_mcse,
                   alpha, beta, gamma, rates, N):
    """Oracle vs measured plus paired difference for each rate and method."""
    reps = estimates.shape[0]
    rows = []
    paired_rows = []
    for k, rate in enumerate(rates):
        for m, method in enumerate(METHODS):
            pair = estimates[:, k, :, m, 0]
            v_o, v_m = pair[:, 0], pair[:, 1]
            so = _stats(v_o, estimates[:, k, 0, m, 1], tau)
            sm = _stats(v_m, estimates[:, k, 1, m, 1], tau)
            common = np.isfinite(v_o) & np.isfinite(v_m)
            d = v_m[common] - v_o[common]
            if len(d) == 0:
                raise ValueError('No valid oracle/measured replicate pairs')
            for label, stat in (('oracle', so), ('measured', sm)):
                rows.append(dict(alpha=alpha, beta=beta, b_degree=gamma,
                                 edge_deletion=float(rate), estimator=method,
                                 weights=label, replicates=reps,
                                 n_villages=N, truth=tau, truth_mcse=tau_mcse,
                                 zero_pair_probability_rate=float(zero[:, k].sum()/(reps*N)),
                                 newly_zero_pair_probability_rate=float(zero_positive[:, k].sum()/(reps*N)),
                                 **stat))
            paired_rows.append(dict(
                alpha=alpha, beta=beta, b_degree=gamma, edge_deletion=float(rate),
                estimator=method, replicates=reps, n_valid_pairs=len(d),
                incremental_bias=float(np.mean(d)),
                incremental_bias_mcse=float(np.std(d, ddof=1)/np.sqrt(len(d))) if len(d)>1 else np.nan,
                mean_absolute_paired_shift=float(np.mean(np.abs(d))),
                paired_shift_sd=float(np.std(d, ddof=1)) if len(d)>1 else np.nan,
                paired_shift_rmse=float(np.sqrt(np.mean(d**2))),
                oracle_variance=so['variance'], measured_variance=sm['variance'],
                variance_change=sm['variance']-so['variance'],
                oracle_bias=so['bias'], measured_bias=sm['bias'],
                oracle_rmse=so['rmse'], measured_rmse=sm['rmse'],
                rmse_change=sm['rmse']-so['rmse'],
                oracle_coverage_90=so['coverage_90'], measured_coverage_90=sm['coverage_90'],
                zero_pair_probability_rate=float(zero[:, k].sum()/(reps*N)),
                newly_zero_pair_probability_rate=float(zero_positive[:, k].sum()/(reps*N)),
            ))
    return pd.DataFrame(rows), pd.DataFrame(paired_rows)


def plot_results(summary, paired, output):
    """Plot pooled descriptive curves; per-cell CSVs retain full comparisons."""
    # Macro-averages across simulation cells are descriptive, not pooled estimates.
    g = paired.groupby(['edge_deletion', 'estimator'], as_index=False).agg(
        incremental_bias=('incremental_bias', 'mean'),
        variance_change=('variance_change', 'mean'),
        rmse_change=('rmse_change', 'mean'),
        newly_zero_pair_probability_rate=('newly_zero_pair_probability_rate', 'mean'),
    )
    g.to_csv(output/'extension_aggregate_paired.csv', index=False)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for name in ('HT', 'Hajek'):
        d = g[g.estimator==name].sort_values('edge_deletion')
        for ax, col in zip(axes, ('incremental_bias','variance_change','rmse_change')):
            ax.plot(d.edge_deletion*100, d[col], 'o-', label=name)
    for ax, title in zip(axes, ('Incremental bias (measured - oracle)',
                                'Variance change (measured - oracle)',
                                'RMSE change (measured - oracle)')):
        ax.axhline(0, c='gray', lw=0.8)
        ax.set_title(title)
        ax.set_xlabel('Directed edge deletion (%)')
        ax.grid(alpha=0.2)
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(output/'extension_paired.png', dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for method in ('HT', 'Hajek'):
        for weights, linestyle in (('oracle', '--'), ('measured', '-')):
            d = summary[(summary.estimator==method)&(summary.weights==weights)]
            d = d.groupby('edge_deletion', as_index=False)[['rmse','coverage_90']].mean().sort_values('edge_deletion')
            axes[0].plot(d.edge_deletion*100, d.rmse, linestyle, marker='o', label=f'{method} {weights}')
            axes[1].plot(d.edge_deletion*100, d.coverage_90, linestyle, marker='o', label=f'{method} {weights}')
    axes[0].set_title('RMSE relative to independent MC truth')
    axes[1].set_title('Nominal 90% CI coverage')
    axes[1].axhline(.9, ls=':', color='gray')
    for ax in axes:
        ax.set_xlabel('Directed edge deletion (%)')
        ax.grid(alpha=0.2)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output/'extension_absolute.png', dpi=160)
    plt.close(fig)


def pack_villages(villages):
    M = len(villages)
    maxn = max(v[1].shape[0] for v in villages)
    A = np.zeros((M,maxn,maxn), np.uint8)
    Din = np.zeros((M,maxn), np.int64)
    Dout = np.zeros((M,maxn), np.int64)
    Q = np.zeros((M,maxn), np.float64)
    Norm = np.zeros(M, np.float64)
    Ns = np.zeros(M, np.int64)
    for i, (_,a,q,di,do,norm) in enumerate(villages):
        n = a.shape[0]
        A[i,:n,:n] = a
        Din[i,:n] = di
        Dout[i,:n] = do
        Q[i,:n] = q
        Norm[i] = norm
        Ns[i] = n
    return A,Din,Dout,Q,Norm,Ns


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--raw', type=Path, default=Path('data/0422allinforawnet.dta'))
    p.add_argument('--replicates', type=int, default=5000, help='Main independent simulation replicates per cell')
    p.add_argument('--truth-replicates', type=int, default=20000, help='Independent reference MC replicates per cell')
    p.add_argument('--n-villages', type=int, default=150)
    p.add_argument('--seed', type=int, default=20261102)
    p.add_argument('--degree-effect', type=float, nargs='+', default=list(DEGREE_EFFECTS))
    p.add_argument('--rates', type=float, nargs='+', default=list(RATES))
    p.add_argument('--alphas', type=float, nargs='+', default=list(ALPHAS))
    p.add_argument('--betas', type=float, nargs='+', default=list(BETAS))
    p.add_argument('--output', type=Path, default=Path('outputs/extension_revised'))
    p.add_argument('--resume', action='store_true')
    p.add_argument('--save-replicates', action='store_true')
    args=p.parse_args()
    if args.replicates<2 or args.truth_replicates<2:
        p.error('Both --replicates and --truth-replicates must be at least 2')
    if any(rate<0 or rate>1 for rate in args.rates):
        p.error('All deletion rates must be between 0 and 1')
    if not args.rates or not args.alphas or not args.betas or not args.degree_effect:
        p.error('Rates and parameter grids cannot be empty')
    args.output.mkdir(parents=True,exist_ok=True)
    villages,diag=load_author_villages(args.raw)
    diag.to_csv(args.output/'network_diagnostics.csv',index=False)
    A,Din,Dout,Q,Norm,Ns=pack_villages(villages)
    if args.n_villages<1 or args.n_villages>len(Ns):
        p.error('n-villages must be between 1 and the number of available villages')
    # Resume only with identical configuration; avoids combining incompatible output files.
    config={
        'raw_path': str(args.raw.resolve()), 'raw_size': args.raw.stat().st_size,
        'replicates': args.replicates, 'truth_replicates': args.truth_replicates,
        'n_villages': args.n_villages, 'seed': args.seed,
        'degree_effect': list(args.degree_effect), 'rates': list(args.rates),
        'alphas': list(args.alphas), 'betas': list(args.betas),
    }
    config_file=args.output/'run_config.json'
    if args.resume and config_file.exists():
        if json.loads(config_file.read_text()) != config:
            raise ValueError('Resume config differs. Choose a new output directory.')
    elif config_file.exists():
        raise FileExistsError('Output directory already contains a run. Use --resume or choose a new --output.')
    config_file.write_text(json.dumps(config,indent=2))
    cells=[(float(g),float(a),float(b)) for g in args.degree_effect for a in args.alphas for b in args.betas]
    summary_files=[]
    paired_files=[]
    for ix,(gamma,alpha,beta) in enumerate(cells):
        stem=f'cell_degree_{gamma:g}_alpha_{alpha:g}_beta_{beta:g}'
        fs=args.output/f'{stem}_summary.csv'
        fp=args.output/f'{stem}_paired.csv'
        ft=args.output/f'{stem}_truth.csv'
        summary_files.append(fs)
        paired_files.append(fp)
        if args.resume and fs.is_file() and fp.is_file() and ft.is_file() and (not args.save_replicates or (args.output/f'{stem}_replicates.npz').exists()):
            print(f'[{ix+1}/{len(cells)}] skip {stem}',flush=True)
            continue
        # Separate RNG streams for estimator and independent reference runs.
        truth_d, truth_nom, truth_rand=independent_reference(
            A,Din,Dout,Q,Ns,alpha,gamma,beta,args.n_villages,
            args.truth_replicates,args.seed+1000003+ix*1009,
        )
        tau=float(truth_d.mean())
        tau_mcse=float(truth_d.std(ddof=1)/np.sqrt(args.truth_replicates))
        pd.DataFrame([dict(alpha=alpha,beta=beta,b_degree=gamma,
                          truth=tau,truth_mcse=tau_mcse,truth_replicates=args.truth_replicates,
                          onehop=float(truth_nom.mean()),random=float(truth_rand.mean()))]).to_csv(ft,index=False)
        est,zeros,newzeros=run_cell(A,Din,Dout,Q,Norm,Ns,alpha,gamma,beta,
                                    args.n_villages,args.replicates,
                                    np.asarray(args.rates,dtype=np.float64), args.seed+ix*1009)
        s,d=summarize_cell(est,zeros,newzeros,tau,tau_mcse,alpha,beta,gamma,args.rates,args.n_villages)
        s.to_csv(fs,index=False)
        d.to_csv(fp,index=False)
        if args.save_replicates:
            np.savez_compressed(args.output/f'{stem}_replicates.npz',
                 estimate=est[:,:,:,:,0], nominal_se=est[:,:,:,:,1],
                 zero_probability_counts=zeros, newly_zero_probability_counts=newzeros,
                 rates=np.asarray(args.rates), estimators=np.asarray(METHODS),
                 weights=np.asarray(WEIGHT_TYPES), truth=tau,truth_mcse=tau_mcse)
        print(f'[{ix+1}/{len(cells)}] done degree={gamma:g}, alpha={alpha:g}, beta={beta:g}; truth={tau:.5f} +/- {tau_mcse:.5f} MCSE',flush=True)
    summary=pd.concat([pd.read_csv(f) for f in summary_files],ignore_index=True)
    paired=pd.concat([pd.read_csv(f) for f in paired_files],ignore_index=True)
    summary.to_csv(args.output/'summary.csv',index=False)
    paired.to_csv(args.output/'paired_comparison.csv',index=False)
    plot_results(summary,paired,args.output)
    print('FINISHED',args.output.resolve(),flush=True)


if __name__=='__main__':
    main()
