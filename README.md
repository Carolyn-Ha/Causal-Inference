# Causal-Inference

Course project repository for **S&DS 6165: Topics in Causal Inference**.

## Project 1: Optimal Trimming and Overlap

This project studies sample trimming based on propensity scores, following Crump et al. (2009), *Dealing with Limited Overlap in Estimation of Average Treatment Effects*.

We focus on two tasks:

### Replication

We replicate selected theoretical and empirical results from Crump et al. (2009), including the optimal propensity-score trimming rule and its effect on estimation variance.

- `Replication1/`: Replication of Table 1 and related simulations
- `Replication2/`: Replication of additional empirical results

### Extension

We extend the original framework to **heteroskedastic outcomes**, where outcome variance can depend on the propensity score.

The extension compares:
- Standard propensity-score trimming
- Optimal symmetric trimming
- Variance-aware trimming

See `Extension1/` for the simulation code, results, and visualizations.

## Reference

Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009).  
*Dealing with limited overlap in estimation of average treatment effects.*  
Biometrika, 96(1), 187–199.

