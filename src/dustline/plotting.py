"""Diagnostic figures for a sightline result.

- plot_law:  the measured R_V with the per-star data points (R_V vs A_V,
  coloured by distance) and the R_V histogram.
- plot_run:  extinction vs distance: the plx S/N > 5 Gaia stars with their
  per-star A_X from the XP x NewEra fits coloured by fitted T_eff, the running
  median with the 16-84 % band, and the dashed bridge to the red-clump bulge
  column when present.
- plot_column: column mode (high latitude): per-star A_V vs D with the stars
  carrying DESI priors marked, the foreground column with its error, an
  optional reference value (SFD/Planck), and the field map of the column.

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


def plot_column(fit: pd.DataFrame, law: dict, path: str, reference_av: float | None = None,
                title: str = "", min_snr: float = 5.0) -> str:
    """Column-mode figure: A_V vs D (left) and the median-A_V field map (right)."""
    from .ensemble import good_sample

    plt = _plt()
    col = law.get("column") or {}
    s = fit[good_sample(fit, min_snr=min_snr)]
    spec = s["spec_prior"].astype(bool) if "spec_prior" in s else pd.Series(False, index=s.index)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw=dict(width_ratios=[1.7, 1]))

    a = ax[0]
    a.errorbar(s.D_kpc, s.av, xerr=s.D_err, yerr=s.av_err, fmt="none", ecolor="0.85", lw=0.5,
               zorder=0)
    sc = a.scatter(s.D_kpc[~spec], s.av[~spec], c=s.teff[~spec], s=14, cmap="RdYlBu",
                   vmin=3500, vmax=7500, marker="o", zorder=2,
                   label=f"XP + photometry only (N = {int((~spec).sum())})")
    if spec.any():
        a.scatter(s.D_kpc[spec], s.av[spec], c=s.teff[spec], s=26, cmap="RdYlBu",
                  vmin=3500, vmax=7500, marker="s", edgecolors="k", linewidths=0.5, zorder=3,
                  label=f"with DESI T$_{{eff}}$/log g/[Fe/H] prior (N = {int(spec.sum())})")
    if col.get("n"):
        d0 = col["d_min_kpc"]
        h = col.get("clean", col)
        a.axvline(d0, color="0.5", ls=":", lw=1)
        a.axhline(h["av"], color="k", lw=1.8, zorder=4,
                  label=f"column A$_V$ = {h['av']:.3f} $\\pm$ {h['av_err']:.3f} "
                        f"(T$_{{eff}}$ $\\geq$ {col.get('teff_clean', 5500):g} K, "
                        f"D > {d0:g} kpc, N = {h['n']})")
        a.axhspan(h["av"] - h["av_mad"], h["av"] + h["av_mad"], color="k", alpha=0.08,
                  label=f"star-to-star MAD {h['av_mad']:.3f}")
        if col.get("teff_cool", {}).get("n"):
            a.axhline(col["teff_cool"]["av"], color="0.4", lw=1, ls="-.", zorder=4,
                      label=f"cool stars (< {col.get('teff_clean', 5500):g} K): "
                            f"{col['teff_cool']['av']:.3f} (template systematic)")
    if reference_av is not None:
        a.axhline(reference_av, color="tab:red", ls="--", lw=1.5, zorder=4,
                  label=f"reference A$_V$ = {reference_av:.3f}")
    a.set_xscale("log")
    a.set_xlim(max(0.05, s.D_kpc.min() * 0.8), s.D_kpc.max() * 1.3)
    a.set_ylim(-0.05, max(0.4, float(np.nanpercentile(s.av, 98)) * 1.3))
    a.set_xlabel("D (kpc, 1/parallax)")
    a.set_ylabel("A$_V$ (per star)")
    a.grid(alpha=0.3, which="both")
    a.legend(fontsize=8, loc="upper left")
    plt.colorbar(sc, ax=a, label="T$_{eff}$ (K, fit)")

    b = ax[1]
    m = col.get("map")
    if m and m.get("cells"):
        n = m["ncell"]
        grid = np.full((n, n), np.nan)
        for k, c in enumerate(m["cells"]):
            if c["av"] is not None:
                grid[k // n, k % n] = c["av"]
        hw = m["half_width_deg"] * 60.0
        im = b.imshow(grid, origin="lower", extent=(-hw, hw, -hw, hw), cmap="viridis")
        for k, c in enumerate(m["cells"]):
            if c["n"]:
                b.text(c["x"] * 60, c["y"] * 60, f"{c['n']}", ha="center", va="center",
                       fontsize=7, color="w")
        b.set_xlabel("$\\Delta$RA cos Dec (arcmin)")
        b.set_ylabel("$\\Delta$Dec (arcmin)")
        b.set_title("median A$_V$ per cell (N stars)", fontsize=10)
        plt.colorbar(im, ax=b, label="A$_V$")
    else:
        b.axis("off")
    fig.suptitle(title or "foreground column", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path
