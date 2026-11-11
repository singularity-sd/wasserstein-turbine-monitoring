# Wasserstein distances for wind-turbine power-curve monitoring: a pre-specified simulation comparison

Wasserstein distances for wind-turbine power-curve monitoring: a pre-specified simulation comparison, from a Singularity v3 forecast.

**Paper:** [paper.pdf](paper.pdf) · **Page:** [singularitysd.com/research/wasserstein.html](https://singularitysd.com/research/wasserstein.html) · Preprint, not yet peer reviewed.

## Abstract

Wind-turbine performance monitoring compares recent SCADA records with a healthy reference, usually through a binned power curve or a residual statistic. These tools compare means; many faults instead change the *shape* of the distribution of power at a given wind speed. We compare five detectors on simulated 10-minute SCADA data with a fixed false-alarm rate of 1% per one-week window: the mean residual, a binned power-curve chi-square statistic, the Kolmogorov–Smirnov statistic on residuals, the 1-Wasserstein distance between residual distributions, and a 1-Wasserstein distance conditional on wind speed. Four fault families at four severities were specified before the single run: uniform derating, a shift of the power curve, a power cap and intermittent half-power operation. Averaged over the twelve non-trivial scenarios, detection rates were 73% (conditional Wasserstein), 71% (residual Wasserstein), 67% (binned chi-square), 54% (mean residual) and 28% (Kolmogorov–Smirnov). The Wasserstein detectors were best for derating and intermittent faults — for example 81–83% against 46% at a 2% intermittent share — while the binned chi-square was clearly best for a power cap (94% against 62% at a 5% cap). No detector dominates; the conditional Wasserstein distance is the best single choice in this benchmark and costs no more than a binned power curve. Code and results are released.

## How it was made

This study started from the Singularity forecast **T3v3-2026-0024**: a method predicted to enter a new field. It was then tested once, with the plan fixed before looking at the results.

## What is in this repository

| file | what it is |
|---|---|
| `experiment.py` | reproduces the study end to end |
| `fig_example.png` | figure |
| `fig_power.png` | figure |
| `main.tex` | LaTeX source of the paper |
| `paper.pdf` | the paper |
| `results.json` | every number reported in the paper |

## Reproduce

```bash
pip install matplotlib numpy scipy
python experiment.py
```

It writes `results.json` and the figures next to the script.

## Cite

```bibtex
@misc{singularity_wasserstein_turbine_monitoring_2026,
  title  = {Wasserstein distances for wind-turbine power-curve monitoring: a pre-specified simulation comparison},
  author = {Singularity},
  year   = {2026},
  url    = {https://github.com/singularity-sd/wasserstein-turbine-monitoring}
}
```

## License

Code: MIT (`LICENSE`). Paper, text and figures: CC BY 4.0 (`LICENSE-paper.md`). Third-party data keep their own licence where stated.

---

Part of [Singularity](https://singularitysd.com): forecasting where science goes next, then testing it.
