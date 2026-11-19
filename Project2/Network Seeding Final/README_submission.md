# Causal Inference Project — Code and Results

This folder contains the R code, input network data, saved simulation results, and figures accompanying the final project report. The project replicates Figure 4 from Chin et al. (2022), *Evaluating Stochastic Seeding Strategies in Networks*, and studies how network structure affects the performance of one-hop versus random seeding, including preferential-attachment networks and degree-preserving rewiring.

## Folder layout

- `01_*.r`–`06_*.r` and `12_optimized_simulation.r`: core simulation, estimation, and Cai-network loading functions.
- `15_full_simulation_150.r`: Figure 4 Monte Carlo runner (computationally intensive).
- `14_figure4_final_analysis.r`: Figure 4 summary and plots.
- `26_coverage_filter_check.r`: coverage-related sensitivity analysis.
- `16_extension_networks.r`–`21_extension_mechanism_analysis.r`: preferential-attachment network generation, main extension simulations, summaries, robustness, and mechanism diagnostics.
- `22_extension_presentation_figures.r`: illustrative network plots and design table.
- `23_*.r`–`25_*.r`: degree-preserving rewiring and clustering extension.
- `cai_data/`: original network input (`0422allinforawnet.dta`), where redistribution is permitted.
- `results/`: saved simulation results and network pools, in the same relative paths used by the scripts.
- `figures/`: selected final figures accompanying the report.

## Software and working directory

Run the scripts in R with the working directory set to **this folder** (the folder containing the `.r` files). The code uses relative paths for `source()`, data loading, and result files. Main packages include `tidyverse`, `igraph`, `haven`, `patchwork`, and `gridExtra`; see individual scripts for dependencies.

## Reproduce results without rerunning the full Monte Carlo simulations

With the included combined result files in place, the following scripts should reconstruct the relevant summaries and plots without rerunning the full Monte Carlo experiments:

```r
source("14_figure4_final_analysis.r")
source("18_extension_analysis.r")
source("20_extension_robustness_analysis.r")
source("26_coverage_filter_check.r")
```

These scripts depend on the saved `.rds` files under `results/`. Files 18 and 20 should be the submission-ready versions that first load combined results and only fall back to individual checkpoints when needed. The exact outputs and other required input paths are documented within the scripts. Some additional analyses (for example, mechanism diagnostics and the rewiring simulation) perform new simulations when sourced; do not run them merely to inspect existing plots.

## Reproduce the simulations from the beginning

The simulation-generation scripts are separated from the result-analysis scripts. The main stages are:

1. `15_full_simulation_150.r` — replicate Figure 4 using Cai networks (requires `cai_data/`).
2. `16_extension_networks.r` then `17_extension_main_simulation.r` — generate synthetic preferential-attachment network pools and run the main extension.
3. `19_extension_robustness.r` — generate additional independent network pools and run robustness simulations.
4. `21_extension_mechanism_analysis.r` — run targeted weight/outcome diagnostics.
5. `23_*.r`, `24_*.r`, then `25_*.r` — generate and analyze the rewiring extension.

These runs may be computationally expensive. Several scripts use checkpoint files to resume interrupted work. Saved result files are provided so that reviewing the report's findings does not require repeating the complete Monte Carlo experiments.

## Notes

- The selected figures are included for direct comparison with the report.
- The raw Cai dataset is supplied only if permitted by its source's distribution terms.
- The `22_extension_presentation_figures.r` script produces illustrative graphs, not necessarily the exact network realizations used in the main simulations.
- Original development/pilot scripts and most intermediate checkpoints have been excluded from this submission to keep it manageable.
