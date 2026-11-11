"""Wasserstein distance for wind-turbine power-curve degradation: a pre-specified simulation study.

From the Singularity v3 forecast T3v3-2026-0024 (wasserstein distance -> Energy Engineering and Power Technology).
Everything below (turbine, wind, noise, faults, detectors, window lengths, false-alarm rate, number of runs, seed)
was fixed before the first run; the script was run once and its output is reported as is.

Run: python experiment.py  ->  results.json, fig_power.png, fig_example.png
"""

import json
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
RNG = np.random.default_rng(20261105)
PR = 2000.0                    # rated power, kW
V_IN, V_OUT, K, V50 = 3.0, 25.0, 0.9, 9.5
N_REF, N_TEST = 4 * 1008, 1008  # four weeks and one week of 10-minute records
BINS = np.arange(3.0, 25.5, 0.5)
ALPHA = 0.01                   # false-alarm rate per test window
N_NULL, N_ALT = 1000, 300
FAULTS = {
    "derating": [0.01, 0.02, 0.03, 0.05],          # P -> (1 - s) P, e.g. blade erosion
    "curve_shift": [0.1, 0.2, 0.3, 0.5],          # v50 -> v50 + s (m/s), e.g. pitch or yaw misalignment
    "power_cap": [0.02, 0.05, 0.08, 0.12],        # P -> min(P, (1 - s) PR), e.g. a controller limit
    "intermittent": [0.01, 0.02, 0.04, 0.08],     # a fraction s of records at half power, e.g. sporadic faults
}


def power_curve(v, v50=V50):
    p = PR / (1.0 + np.exp(-K * (v - v50)))
    return np.where((v >= V_IN) & (v <= V_OUT), p, 0.0)


def sample(n, fault=None, s=0.0):
    v = stats.weibull_min.rvs(2.0, scale=8.0, size=n, random_state=RNG)
    p = power_curve(v, V50 + (s if fault == "curve_shift" else 0.0))
    if fault == "derating":
        p = p * (1 - s)
    if fault == "power_cap":
        p = np.minimum(p, (1 - s) * PR)
    if fault == "intermittent":
        p = np.where(RNG.random(n) < s, 0.5 * p, p)
    p = p * (1 + RNG.normal(0, 0.04, n)) + RNG.normal(0, 0.01 * PR, n)
    return v, np.clip(p, 0, PR)


class Reference:
    def __init__(self, v, p):
        self.idx = np.digitize(v, BINS)
        self.med = np.array([np.median(p[self.idx == b]) if (self.idx == b).sum() >= 5 else np.nan for b in range(len(BINS) + 1)])
        self.var = np.array([np.var(p[self.idx == b]) if (self.idx == b).sum() >= 5 else np.nan for b in range(len(BINS) + 1)])
        self.byb = {b: np.sort(p[self.idx == b]) for b in range(len(BINS) + 1) if (self.idx == b).sum() >= 5}
        self.res = p - self.med[self.idx]
        ok = np.isfinite(self.res)
        self.res = self.res[ok]


def detectors(ref, v, p):
    idx = np.digitize(v, BINS)
    r = p - ref.med[idx]
    ok = np.isfinite(r)
    r, idx_ok, p_ok = r[ok], idx[ok], p[ok]
    out = {}
    out["mean_residual"] = abs(r.mean()) / (ref.res.std() / np.sqrt(len(r)))
    chi = 0.0
    for b in np.unique(idx_ok):
        m = idx_ok == b
        if m.sum() >= 3 and np.isfinite(ref.var[b]) and ref.var[b] > 0:
            chi += m.sum() * (r[m].mean()) ** 2 / ref.var[b]
    out["binned_chi2"] = chi
    out["ks_residual"] = stats.ks_2samp(r, ref.res).statistic
    out["w1_residual"] = stats.wasserstein_distance(r, ref.res)
    w, tot = 0.0, 0
    for b in np.unique(idx_ok):
        m = idx_ok == b
        if m.sum() >= 3 and b in ref.byb:
            w += m.sum() * stats.wasserstein_distance(p_ok[m], ref.byb[b])
            tot += m.sum()
    out["w1_conditional"] = w / tot
    return out


def main():
    ref = Reference(*sample(N_REF))
    null = {k: [] for k in ["mean_residual", "binned_chi2", "ks_residual", "w1_residual", "w1_conditional"]}
    for _ in range(N_NULL):
        for k, x in detectors(ref, *sample(N_TEST)).items():
            null[k].append(x)
    thr = {k: float(np.quantile(v, 1 - ALPHA)) for k, v in null.items()}
    res = {"protocol": {"N_REF": N_REF, "N_TEST": N_TEST, "ALPHA": ALPHA, "N_NULL": N_NULL, "N_ALT": N_ALT, "faults": FAULTS},
           "thresholds": thr, "power": {}}
    for f, sev in FAULTS.items():
        res["power"][f] = {}
        for s in sev:
            hits = {k: 0 for k in thr}
            for _ in range(N_ALT):
                d = detectors(ref, *sample(N_TEST, f, s))
                for k in thr:
                    hits[k] += d[k] > thr[k]
            res["power"][f][str(s)] = {k: hits[k] / N_ALT for k in thr}
            print(f, s, res["power"][f][str(s)], flush=True)
    (HERE / "results.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    figures(res)


def figures(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = {"mean_residual": "mean residual", "binned_chi2": "binned power curve (chi-square)", "ks_residual": "Kolmogorov-Smirnov, residuals",
             "w1_residual": "Wasserstein, residuals", "w1_conditional": "Wasserstein, conditional on wind speed"}
    titles = {"derating": "Derating (P × (1−s))", "curve_shift": "Curve shift (v50 + s m/s)", "power_cap": "Power cap ((1−s) × rated)",
              "intermittent": "Intermittent half power (share s)"}
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.6))
    for ax, (f, d) in zip(axs, res["power"].items()):
        xs = [float(s) for s in d]
        for k, lab in names.items():
            ax.plot(xs, [d[str(s) if str(s) in d else s][k] for s in d], marker="o", label=lab,
                    lw=2.4 if k == "w1_conditional" else 1.4)
        ax.set_title(titles[f], fontsize=10)
        ax.set_xlabel("severity s")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(alpha=0.3)
    axs[0].set_ylabel(f"detection rate (false alarms {ALPHA:.0%})")
    axs[-1].legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(HERE / "fig_power.png", dpi=160)
    v, p = sample(3000)
    vf, pf = sample(3000, "power_cap", 0.08)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(v, p, s=2, alpha=0.3, label="healthy")
    ax.scatter(vf, pf, s=2, alpha=0.3, label="power cap s = 0.08")
    ax.set_xlabel("wind speed (m/s)")
    ax.set_ylabel("power (kW)")
    ax.legend(markerscale=6)
    fig.tight_layout()
    fig.savefig(HERE / "fig_example.png", dpi=160)


if __name__ == "__main__":
    main()
