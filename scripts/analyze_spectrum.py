"""Analyze the JWST NIRCam transmission spectrum of HD 189733 b.

Data source: Zenodo record 10.5281/zenodo.11459715, "Products for
'Hydrogen sulfide and metal-enriched atmosphere for a Jupiter-mass
exoplanet'", file Data/transit spectrum/HD189733b_NIRCam_spectrum.txt.
Retrieved directly from Zenodo; reproduced unmodified in data/.

This script computes the weighted mean transit depth and compares the
mean depth in two wavelength windows against a nearby continuum window:
one spanning the CO2 band near 4.3-4.6 micron, and one spanning
2.6-3.0 micron where H2O and H2S both contribute. Both statistics are
band-contrast signal-to-noise ratios on a simple two-window comparison,
not molecular detection significances -- a broadband window like this
can contain contributions from more than one absorber, and the
continuum itself is a modeling choice, not a measured baseline. Fu et
al. (2024) derive the actual per-molecule detection significances from
a full retrieval/model comparison: H2O at 13.4 sigma, CO2 at 11.2
sigma, CO at 5 sigma, and H2S at 4.5 sigma. This script's own numbers
are reported alongside those, not as a substitute for them.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import scienceplots  # noqa: F401 (registers 'science' style)
import numpy as np

plt.style.use(["science", "no-latex"])

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"

CO2_BAND = (4.3, 4.6)
H2O_H2S_BAND = (2.6, 3.0)
CONTINUUM_BAND = (3.6, 3.9)

# Fu et al. (2024) full retrieval/model-comparison significances -- a
# different, more rigorous estimator than this script's band contrast.
PAPER_SIGMA = {"H2O": 13.4, "CO2": 11.2, "CO": 5.0, "H2S": 4.5}


def load_spectrum(path: Path):
    wave, depth, err = [], [], []
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 4:
                continue
            w, _werr, d, derr = map(float, parts)
            wave.append(w)
            depth.append(d)
            err.append(derr)
    return np.array(wave), np.array(depth), np.array(err)


def weighted_mean(values: np.ndarray, errors: np.ndarray) -> tuple[float, float]:
    weights = 1.0 / errors**2
    mean = np.sum(values * weights) / np.sum(weights)
    mean_error = np.sqrt(1.0 / np.sum(weights))
    return mean, mean_error


def band_stats(wave, depth, err, band):
    mask = (wave >= band[0]) & (wave <= band[1])
    return weighted_mean(depth[mask], err[mask])


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    wave, depth, err = load_spectrum(DATA_DIR / "nircam_transmission_spectrum.txt")
    order = np.argsort(wave)
    wave, depth, err = wave[order], depth[order], err[order]

    mean_depth, mean_depth_err = weighted_mean(depth, err)
    co2_mean, co2_err = band_stats(wave, depth, err, CO2_BAND)
    h2o_mean, h2o_err = band_stats(wave, depth, err, H2O_H2S_BAND)
    cont_mean, cont_err = band_stats(wave, depth, err, CONTINUUM_BAND)

    co2_excess_ppm = (co2_mean - cont_mean) * 1e6
    co2_sigma = abs(co2_mean - cont_mean) / np.sqrt(co2_err**2 + cont_err**2)
    h2o_excess_ppm = (h2o_mean - cont_mean) * 1e6
    h2o_sigma = abs(h2o_mean - cont_mean) / np.sqrt(h2o_err**2 + cont_err**2)

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        writer.writerow(["n_wavelength_bins", len(wave), "count"])
        writer.writerow(["wavelength_min", f"{wave.min():.3f}", "micron"])
        writer.writerow(["wavelength_max", f"{wave.max():.3f}", "micron"])
        writer.writerow(["weighted_mean_depth", f"{mean_depth*1e6:.1f}", "ppm"])
        writer.writerow(["continuum_mean_depth", f"{cont_mean*1e6:.1f}", "ppm (3.6-3.9 um)"])
        writer.writerow(["co2_band_mean_depth", f"{co2_mean*1e6:.1f}", "ppm (4.3-4.6 um)"])
        writer.writerow(["co2_band_excess", f"{co2_excess_ppm:.1f}", "ppm"])
        writer.writerow(["co2_band_contrast_snr_this_script", f"{co2_sigma:.1f}", "sigma (two-window band contrast, not a molecular detection significance)"])
        writer.writerow(["h2o_h2s_band_mean_depth", f"{h2o_mean*1e6:.1f}", "ppm (2.6-3.0 um)"])
        writer.writerow(["h2o_h2s_band_excess", f"{h2o_excess_ppm:.1f}", "ppm"])
        writer.writerow(["h2o_h2s_band_contrast_snr_this_script", f"{h2o_sigma:.1f}", "sigma (two-window band contrast, not a molecular detection significance)"])
        for molecule, sigma in PAPER_SIGMA.items():
            writer.writerow([f"paper_retrieval_significance_{molecule}", f"{sigma}", "sigma (Fu et al. 2024, full retrieval)"])

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    ax.errorbar(wave, depth * 1e6, yerr=err * 1e6, fmt="o", ms=3, color="#2c5f8a", ecolor="#9fbfd8", elinewidth=0.8, label="HD 189733 b, JWST NIRCam")
    ax.axhline(mean_depth * 1e6, color="#555555", lw=1, ls="--", label="weighted mean")
    ax.axvspan(*H2O_H2S_BAND, color="#1f6f5c", alpha=0.12, label="H2O/H2S-influenced band (2.6-3.0 um)")
    ax.axvspan(*CO2_BAND, color="#c0562a", alpha=0.12, label="CO2 band (4.3-4.6 um)")
    ax.axvspan(*CONTINUUM_BAND, color="#999999", alpha=0.12, label="continuum window (3.6-3.9 um)")
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Transit depth (Rp/Rs)^2 [ppm]")
    ax.set_title("HD 189733 b transmission spectrum (JWST NIRCam)")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "hd189733b_transmission_spectrum.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {FIG_DIR / 'hd189733b_transmission_spectrum.png'}")
    print(f"n={len(wave)}, weighted mean depth = {mean_depth*1e6:.1f} +/- {mean_depth_err*1e6:.2f} ppm")
    print(f"CO2-band contrast (this script) = {co2_excess_ppm:.1f} ppm ({co2_sigma:.1f} sigma band S/N)")
    print(f"H2O/H2S-band contrast (this script) = {h2o_excess_ppm:.1f} ppm ({h2o_sigma:.1f} sigma band S/N)")
    print(f"Fu et al. (2024) retrieval significances: H2O {PAPER_SIGMA['H2O']} sigma, CO2 {PAPER_SIGMA['CO2']} sigma, CO {PAPER_SIGMA['CO']} sigma, H2S {PAPER_SIGMA['H2S']} sigma")


if __name__ == "__main__":
    main()
