"""The public API: Sightline and ExtinctionResult.

    import dustline
    res = dustline.Sightline(ra=275.089125, dec=-18.269389).run()
    res.rv                 # {'median': ..., 'mad': ..., 'n_stars': ...}
    res.extinction("I")    # A_I(D) DataFrame
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import ensemble
from .cache import Workspace
from .catalogs import footprint, gaia
from .userphot import UserPhotometry

DEFAULT_RADIUS_ARCMIN = 5.0
DEFAULT_FILTER = "I"       # Cousins I


class ExtinctionResult:
    """The measured law + dust run for one sightline."""

    def __init__(self, ws: Workspace, law: dict, run_av: pd.DataFrame, plan):
        self._ws = ws
        self.law = law
        self._run_av = run_av
        self.plan = plan

    @property
    def stars(self) -> pd.DataFrame:
        """The per-star fit table (lazy-loaded from the workspace)."""
        return pd.read_csv(self._ws.path("xp_stars.csv"))

    @property
    def clump_anchor(self) -> dict | None:
        """The red-clump bulge anchor used for the bridge (None off the bulge)."""
        return self.law.get("clump_anchor")

    @property
    def mode(self) -> str:
        """"law" (measured R_V + run) or "column" (high latitude: foreground column)."""
        return self.law.get("mode", "law")

    @property
    def column(self) -> dict | None:
        """Column mode: the foreground A_V of the field (ensemble.foreground_column)."""
        return self.law.get("column")

    def plots(self, band: str = DEFAULT_FILTER, directory: str | None = None,
              reference_av: float | None = None) -> dict:
        """Write the diagnostic figures; returns {'law' | 'column': path, 'run': path}.

        - law (law mode): per-star R_V vs A_V (the law sample) + R_V histogram;
        - column (column mode): per-star A_V vs D with the foreground column and
          its field map (reference_av: e.g. the SFD/Planck A_V, drawn for comparison);
        - run: A_X vs D for the plx S/N > 5 stars coloured by fitted T_eff, the
          running median with the 16-84 % band, and the dashed red-clump bridge.
        """
        from pathlib import Path

        from . import plotting

        d = Path(directory) if directory else self._ws.dir
        d.mkdir(parents=True, exist_ok=True)
        fit = self.stars
        tag = (f"({self._ws.ra:.4f}, {self._ws.dec:.4f})")
        out = {"run": str(d / f"extinction_run_{band}.png")}
        if self.mode == "column":
            out["column"] = str(d / "column_av.png")
            plotting.plot_column(fit, self.law, out["column"], reference_av=reference_av,
                                 title=f"{tag}: foreground column")
        else:
            out["law"] = str(d / "law_rv.png")
            plotting.plot_law(fit, self.law, out["law"],
                              title=f"{tag}: R$_V$ = {self.law['rv']:.2f} $\\pm$ "
                                    f"{self.law['rv_mad']:.2f} ({self.law['n_stars']} stars)")
        ratio = (self.law["ratios_av"].get(band)
                 or ensemble.band_ratio_at_rv(band, self.law["rv"]))
        plotting.plot_run(fit, self.law, self._run_av, out["run"], band=band, ratio=ratio,
                          title=f"{tag}: extinction vs distance ({band})")
        for k in out.values():
            print(f"wrote {k}")
        return out

    @property
    def rv(self) -> dict:
        """Measured R_V: median, star-to-star MAD, weighted mean, N stars."""
        return {k: self.law[k] for k in
                ("rv", "rv_mad", "rv_mean", "rv_mean_err", "n_stars", "rv_16", "rv_84")}

    def extinction(self, band: str = DEFAULT_FILTER) -> pd.DataFrame:
        """A_X(D): columns D_kpc, A_med, A_16, A_84, n_stars, bridged.

        The run is measured in A_V (per-star, SED-marginalised); the requested
        band is reached through the band-integrated A_X/A_V ratio of a
        reference red-giant SED under G23 at the measured R_V. For bands that
        were IN the fit, the per-star measured ratio is used instead.
        """
        if band in self.law["ratios_av"]:
            ratio = self.law["ratios_av"][band]
        else:
            ratio = ensemble.band_ratio_at_rv(band, self.law["rv"])
        out = self._run_av.copy()
        for c in ("AV_med", "AV_16", "AV_84"):
            out[c.replace("AV", "A")] = out[c] * ratio
        out["band"] = band
        out["ratio_av"] = ratio
        return out[["D_kpc", "A_med", "A_16", "A_84", "n_stars", "bridged", "band", "ratio_av"]]

    def save(self, path: str, band: str = DEFAULT_FILTER) -> None:
        self.extinction(band).round(4).to_csv(path, index=False)

    def __repr__(self) -> str:  # pragma: no cover
        d = self._run_av
        if self.mode == "column" and self.column:
            c = self.column["clean"]
            return (f"ExtinctionResult(column A_V {c['av']:.3f} +/- {c['av_err']:.3f} "
                    f"({c['n']} F/G stars beyond {self.column['d_min_kpc']:g} kpc; "
                    f"R_V {self.law['rv']:g} assumed), "
                    f"run {d.D_kpc.min():.1f}-{d.D_kpc.max():.1f} kpc)")
        return (f"ExtinctionResult(R_V {self.law['rv']:.2f} +/- {self.law['rv_mad']:.2f} "
                f"({self.law['n_stars']} stars), run {d.D_kpc.min():.1f}-{d.D_kpc.max():.1f} kpc)")


class Sightline:
    """A line of sight: give RA/Dec (deg, ICRS), run the pipeline, get the result.

    Parameters
    ----------
    ra, dec : float          position in degrees
    radius_arcmin : float    field radius (default 5'); larger = more stars, slower
    photometry : UserPhotometry | None    your own calibrated photometry
    prefer_deep : bool       use DECaPS/VVV instead of PS1/2MASS where available
    spectro_priors : bool    use DESI MWS T_eff/log g/[Fe/H] as per-star template
                             priors where the footprint has them (high latitude)
    freeze_offsets : bool    hold the photometric zero points at 0 instead of
                             iterating them (the A/B control in column mode)
    """

    def __init__(self, ra: float, dec: float, radius_arcmin: float = DEFAULT_RADIUS_ARCMIN,
                 photometry: UserPhotometry | None = None, prefer_deep: bool = True,
                 spectro_priors: bool = True, freeze_offsets: bool = False):
        self.ra, self.dec = float(ra), float(dec)
        self.radius_arcmin = float(radius_arcmin)
        self.photometry = photometry
        self.plan = footprint.plan(self.ra, self.dec, prefer_deep=prefer_deep)
        if not spectro_priors:
            self.plan.spectro = "none"
        self.freeze_offsets = bool(freeze_offsets)
        config = dict(optical=self.plan.optical, nir=self.plan.nir,
                      user_phot=str(photometry.path) if photometry else None)
        # non-default options only, so the cache keys of plain law-mode runs are unchanged
        if self.plan.spectro != "none":
            config["spectro"] = self.plan.spectro
        if self.plan.mode != "law":
            config["mode"] = self.plan.mode
        if self.freeze_offsets:
            config["freeze_offsets"] = True
        self.ws = Workspace(self.ra, self.dec, self.radius_arcmin, config)

    def run(self, force: bool = False, chunk: int = 200) -> ExtinctionResult:
        """The full pipeline (each stage cached; the XP fetch is resumable).

        First run per sightline: hours (Gaia XP fetch dominates). After: instant.
        """
        from . import bands as bands_mod
        from . import fit as fit_mod

        result_path = self.ws.path("result.json")
        run_path = self.ws.path("dust_run_av.csv")
        if result_path.exists() and run_path.exists() and not force:
            law = json.loads(result_path.read_text())
            return ExtinctionResult(self.ws, law, pd.read_csv(run_path), self.plan)

        print(f"dustline: ({self.ra:.5f}, {self.dec:.5f}), r = {self.radius_arcmin}'")
        for note in self.plan.notes:
            print(f"  - {note}")

        # 1. Gaia cone + XP spectra (fetch skipped when calibrated spectra exist,
        #    e.g. a seeded or completed workspace)
        g = gaia.cone(self.ws, force=force)
        if not self.ws.has("xp_sampled.npz") or force:
            gaia.fetch_xp(self.ws, g, chunk=chunk)
        xp = gaia.calibrate_xp(self.ws, force=force)

        # 2. photometry assembly
        user = self.photometry.matched_to_gaia(g, self.ra, self.dec) if self.photometry else None
        stars, band_list = bands_mod.assemble(self.ws, g, self.plan, user_phot=user, force=force)
        if not band_list:
            raise RuntimeError("no photometric bands available: no survey coverage and no "
                               "user photometry - the XP-only fit is biased and refused")

        # 3. per-star fits with zero-point iteration
        #    (column mode: R_V held at the assumed value - unconstrained at A_V ~ 0.1)
        fit = fit_mod.run_fit_with_offsets(
            self.ws, stars, xp, band_list, force=force,
            spectro_priors=(self.plan.spectro != "none"), freeze_offsets=self.freeze_offsets,
            rv_fixed=ensemble.RV_ASSUMED if self.plan.mode == "column" else None)

        # 4. ensemble products: the measured law, or the assumed one where there
        #    is no reddening to measure it (column mode / too few A_V >= 2 stars)
        if self.plan.mode == "column":
            law = ensemble.default_law(band_list, reason="high-latitude column mode")
        else:
            try:
                law = ensemble.measure_law(fit, band_list)
            except RuntimeError as e:
                print(f"law not measured: {e}")
                law = ensemble.default_law(band_list, reason=str(e).split(" - ")[0])
        run_av = ensemble.dust_run(fit)
        law["bands"] = band_list
        law["n_fitted"] = int(len(fit))
        law["n_spec_prior"] = int(fit["spec_prior"].sum()) if "spec_prior" in fit else 0
        law["survey_plan"] = dict(optical=self.plan.optical, nir=self.plan.nir,
                                  spectro=self.plan.spectro, mode=self.plan.mode)
        off_path = self.ws.path("phot_offsets.json")
        law["phot_offsets"] = json.loads(off_path.read_text()) if off_path.exists() else {}
        law["offsets_frozen"] = self.freeze_offsets
        if self.plan.mode == "column":
            law["column"] = ensemble.foreground_column(fit)

        # 5. red-clump bulge bridge (bulge window only)
        if self.plan.in_bulge_window:
            from . import clump

            nir_prefix = "VISTA" if self.plan.nir == "vvv" else "2MASS"
            anchor = clump.find_clump(stars, law, f"{nir_prefix}_J", f"{nir_prefix}_Ks")
            if anchor:
                law["clump_anchor"] = anchor
                run_av = ensemble.bridge_to_clump(run_av, anchor)
                print(f"red-clump anchor: A_V = {anchor['AV_column']:.2f} +/- "
                      f"{anchor['AV_column_err']:.2f} at D_RC = {anchor['D_RC']:.1f} kpc "
                      f"(E(J-Ks) {anchor['E_JK']:.2f}, {anchor['n_window']} stars) - run bridged")
            else:
                print("bulge window but no credible red-clump anchor - run stops at "
                      "the last parallax bin")

        json.dump(law, open(result_path, "w"), indent=1)
        run_av.round(4).to_csv(run_path, index=False)
        rv = law["rv"]
        if law.get("mode") == "column":
            c = law["column"]
            h = c["clean"]
            print(f"dustline: foreground column A_V = {h['av']:.3f} +/- {h['av_err']:.3f} "
                  f"(MAD {h['av_mad']:.3f}, {h['n']} F/G stars beyond {c['d_min_kpc']:g} kpc; "
                  f"all {c['n']} stars: {c['av']:.3f}; {law['n_spec_prior']} with DESI priors); "
                  f"R_V = {rv:g} assumed")
        else:
            print(f"dustline: R_V = {rv:.2f} +/- {law['rv_mad']:.2f} ({law['n_stars']} stars); "
                  f"run to {run_av.D_kpc.max():.1f} kpc")
        return ExtinctionResult(self.ws, law, run_av, self.plan)


def _to_float(x):
    return float(np.asarray(x))
