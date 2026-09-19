"""Diagnostic figures for a sightline result.

- plot_law:  the measured R_V with the per-star data points (R_V vs A_V,
  coloured by distance) and the R_V histogram.
- plot_run:  extinction vs distance: the plx S/N > 5 Gaia stars with their
  per-star A_X from the XP x NewEra fits coloured by fitted T_eff, the running
  median with the 16-84 % band, and the dashed bridge to the red-clump bulge
  column when present.

matplotlib is an optional dependency (pip install dustline[plot]).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_law(fit: pd.DataFrame, law: dict, path: str, title: str = "") -> str:
    """R_V figure: per-star R_V vs A_V (law sample) + histogram."""
    from .ensemble import law_sample

    plt = _plt()
    w = fit[law_sample(fit, min_av=law.get("min_av", 2.0))]
    rv_med, rv_mad = law["rv"], law["rv_mad"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw=dict(width_ratios=[1.6, 1]))

    sc = ax[0].scatter(w.av, w.rv, c=w.D_kpc, s=14, cmap="viridis", vmin=0.5, vmax=8)
    ax[0].errorbar(w.av, w.rv, yerr=w.rv_err, fmt="none", ecolor="0.75", lw=0.5, zorder=0)
    ax[0].axhline(rv_med, color="k", lw=1.5,
                  label=f"median R$_V$ = {rv_med:.2f} $\\pm$ {rv_mad:.2f} (MAD)")
    ax[0].axhspan(rv_med - rv_mad, rv_med + rv_mad, color="k", alpha=0.10)
    ax[0].set_xlabel("A$_V$ (per star)")
    ax[0].set_ylabel("R$_V$ (G23, per star)")
    ax[0].set_ylim(2.1, 5.8)
    ax[0].legend(fontsize=9, loc="upper right")
    ax[0].grid(alpha=0.3)
    plt.colorbar(sc, ax=ax[0], label="D (kpc)")

    ax[1].hist(w.rv, bins=np.arange(2.2, 5.7, 0.25), color="C0", alpha=0.75)
    ax[1].axvline(rv_med, color="k", lw=1.5)
    ax[1].axvspan(rv_med - rv_mad, rv_med + rv_mad, color="k", alpha=0.10)
    ax[1].set_xlabel("R$_V$")
    ax[1].set_ylabel("stars")
    ax[1].grid(alpha=0.3)

    fig.suptitle(title or f"measured extinction law: R$_V$ = {rv_med:.2f} $\\pm$ {rv_mad:.2f}, "
                          f"{law['n_stars']} stars with A$_V \\geq$ {law.get('min_av', 2.0):g}",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def plot_run(fit: pd.DataFrame, law: dict, run: pd.DataFrame, path: str,
             band: str = "I", ratio: float | None = None, title: str = "",
             min_snr: float = 5.0) -> str:
    """Extinction-distance figure in band X: plx S/N > min_snr stars (A_X per star,
    coloured by fitted T_eff), running median + 16-84 % band, dashed clump bridge."""
    from .ensemble import good_sample

    plt = _plt()
    if ratio is None:
        ratio = law["ratios_av"].get(band) or 1.0
    s = fit[good_sample(fit, min_snr=min_snr)]
    # per-star A_X: the fit's own band extinction when available, else av * ratio
    if f"A_{band}" in s.columns:
        A = s[f"A_{band}"].values
        Aerr = s.get(f"A_{band}_err", pd.Series(np.zeros(len(s)))).values
    else:
        A = s.av.values * ratio
        Aerr = s.av_err.values * ratio

    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    ax.errorbar(s.D_kpc, A, xerr=s.D_err, yerr=Aerr, fmt="none",
                ecolor="0.8", lw=0.5, zorder=0)
    sc = ax.scatter(s.D_kpc, A, c=s.teff, s=14, cmap="RdYlBu", vmin=3500, vmax=7500,
                    zorder=2, label=f"Gaia stars, plx S/N > {min_snr:g} (N = {len(s)})")

    meas = run[~run.bridged]
    ax.plot(meas.D_kpc, meas.AV_med * ratio, "k-", lw=2, zorder=3, label="running median")
    ax.fill_between(meas.D_kpc, meas.AV_16 * ratio, meas.AV_84 * ratio,
                    color="k", alpha=0.12, label="16-84 %")
    br = run[run.bridged]
    if len(br):
        # connect from the last measured point
        Db = np.concatenate([[meas.D_kpc.iloc[-1]], br.D_kpc.values])
        Ab = np.concatenate([[meas.AV_med.iloc[-1]], br.AV_med.values]) * ratio
        ax.plot(Db, Ab, "k--", lw=2, zorder=3, label="bridge to red-clump bulge column")
        ax.fill_between(br.D_kpc, br.AV_16 * ratio, br.AV_84 * ratio,
                        color="k", alpha=0.06)
        anchor = law.get("clump_anchor")
        if anchor:
            ax.errorbar([anchor["D_RC"]], [anchor["AV_column"] * ratio],
                        yerr=[anchor["AV_column_err"] * ratio], fmt="s", color="tab:red",
                        ms=8, mec="k", zorder=4,
                        label=f"red clump: A = {anchor['AV_column'] * ratio:.2f} "
                              f"at {anchor['D_RC']:.1f} kpc")

    ax.set_xscale("log")
    ax.set_xlim(0.3, max(12.0, run.D_kpc.max() * 1.1))
    ax.set_ylim(0, max(np.nanpercentile(A, 99) * 1.3, (run.AV_84.max() * ratio) * 1.15))
    ax.set_xlabel("D (kpc, 1/parallax)")
    label = f"A$_{{{band}}}$" if len(band) <= 2 else f"A ({band})"
    ax.set_ylabel(label)
    ax.legend(fontsize=8.5, loc="upper left")
    ax.grid(alpha=0.3, which="both")
    plt.colorbar(sc, ax=ax, label="T$_{eff}$ (K, fit)")
    ax.set_title(title or f"extinction vs distance ({band}); "
                          f"R$_V$ = {law['rv']:.2f} $\\pm$ {law['rv_mad']:.2f}", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path
