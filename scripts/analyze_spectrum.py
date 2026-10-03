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
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
MODEL_DIR = DATA_DIR / "models"

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


def load_model(path: Path) -> tuple[np.ndarray, np.ndarray]:
    table = np.loadtxt(path)
    table = table[np.all(np.isfinite(table), axis=1)]
    wavelength, depth = table[:, 0], table[:, 1]
    # The grid-family R=300 full model stores micron values multiplied by 1e6.
    if np.median(wavelength) > 1e4:
        wavelength = wavelength / 1e6
    order = np.argsort(wavelength)
    return wavelength[order], depth[order]


def compare_model(
    wave: np.ndarray, depth: np.ndarray, error: np.ndarray, model_path: Path
) -> dict[str, float | str]:
    model_wave, model_depth = load_model(model_path)
    if wave.min() < model_wave.min() or wave.max() > model_wave.max():
        raise ValueError(f"model does not cover data range: {model_path}")
    prediction = np.interp(wave, model_wave, model_depth)
    weights = 1.0 / error**2
    offset = float(np.sum((depth - prediction) * weights) / np.sum(weights))
    residual = (depth - prediction - offset) / error
    chi2 = float(np.sum(residual**2))
    return {
        "model": model_path.stem,
        "chi2": chi2,
        "dof": int(wave.size - 1),
        "reduced_chi2": chi2 / (wave.size - 1),
        "fitted_offset_ppm": offset * 1e6,
    }


def model_family_metrics(
    wave: np.ndarray, depth: np.ndarray, error: np.ndarray, family: str
) -> list[dict[str, float | str]]:
    rows = [compare_model(wave, depth, error, path) for path in sorted((MODEL_DIR / family).glob("*.txt"))]
    full = next(row for row in rows if str(row["model"]).startswith("FULL"))
    for row in rows:
        row["delta_chi2_vs_full_median_curve"] = float(row["chi2"]) - float(full["chi2"])
        row["family"] = family
    return rows


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
    model_rows = model_family_metrics(wave, depth, err, "free") + model_family_metrics(
        wave, depth, err, "grid"
    )
    with (FIG_DIR / "model_comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(model_rows[0]))
        writer.writeheader()
        writer.writerows(model_rows)

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

    payload = {
        "source": {
            "zenodo_record": "10.5281/zenodo.11459715",
            "archive_md5": "7b39f9b6fdfc0d19a919a6d27eebde28",
            "n_spectrum_bins": int(wave.size),
        },
        "supplied_median_curve_comparisons": model_rows,
        "inference_boundary": (
            "One fitted vertical offset and diagonal errors. Delta chi-square between supplied "
            "median curves is descriptive and is not a retrieval Bayes factor or molecule significance."
        ),
        "published_retrieval_significance_sigma": PAPER_SIGMA,
    }
    (FIG_DIR / "analysis_summary.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )

    fig, (ax, ax_r, ax_d) = plt.subplots(3, 1, figsize=(10, 11), constrained_layout=True)
    ax.errorbar(wave, depth * 1e6, yerr=err * 1e6, fmt="o", ms=3, color="#2c5f8a", ecolor="#9fbfd8", elinewidth=0.8, label="HD 189733 b, JWST NIRCam")
    full_predictions = {}
    for family, color in (("free", "#b44b2a"), ("grid", "#2a7f62")):
        row = next(item for item in model_rows if item["family"] == family and str(item["model"]).startswith("FULL"))
        path = next((MODEL_DIR / family).glob("FULL*.txt"))
        model_wave, model_depth = load_model(path)
        prediction = np.interp(wave, model_wave, model_depth) + float(row["fitted_offset_ppm"]) / 1e6
        full_predictions[family] = prediction
        ax.plot(wave, prediction * 1e6, color=color, lw=1.6, label=f"{family} full median + offset")
    ax.axvspan(*H2O_H2S_BAND, color="#1f6f5c", alpha=0.12, label="H2O/H2S-influenced band (2.6-3.0 um)")
    ax.axvspan(*CO2_BAND, color="#c0562a", alpha=0.12, label="CO2 band (4.3-4.6 um)")
    ax.axvspan(*CONTINUUM_BAND, color="#999999", alpha=0.12, label="continuum window (3.6-3.9 um)")
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Transit depth (Rp/Rs)^2 [ppm]")
    ax.set_title("A  JWST/NIRCam spectrum and supplied full-model medians")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax.grid(alpha=0.25)

    for family, color in (("free", "#b44b2a"), ("grid", "#2a7f62")):
        ax_r.plot(wave, (depth - full_predictions[family]) / err, marker="o", ms=2.5, lw=0.8, color=color, label=family)
    ax_r.axhline(0, color="#555", lw=1)
    ax_r.set(title="B  Standardized residuals", ylabel="(data − model) / σ", xlabel="Wavelength [micron]")
    ax_r.legend(frameon=False)
    ax_r.grid(alpha=0.2)

    labels = ["cloud", "−CO", "−CO2", "−H2O", "−H2S"]
    x = np.arange(len(labels))
    for index, (family, color) in enumerate((("free", "#b44b2a"), ("grid", "#2a7f62"))):
        family_rows = {str(row["model"]): row for row in model_rows if row["family"] == family}
        keys = ["Cloud_only_TD_median_result", "xCO_TD_median_result", "xCO2_TD_median_result", "xH2O_TD_median_result", "xH2S_TD_median_result"]
        values = [float(family_rows[key]["delta_chi2_vs_full_median_curve"]) for key in keys]
        ax_d.bar(x + (index - 0.5) * 0.34, values, width=0.34, color=color, label=family)
    ax_d.set_xticks(x, labels)
    ax_d.set(title="C  Shape mismatch of supplied ablation medians", ylabel="Δχ² vs full median (descriptive)")
    ax_d.legend(frameon=False)
    ax_d.grid(axis="y", alpha=0.2)
    fig.savefig(FIG_DIR / "hd189733b_transmission_spectrum.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {FIG_DIR / 'hd189733b_transmission_spectrum.png'}")
    print(f"n={len(wave)}, weighted mean depth = {mean_depth*1e6:.1f} +/- {mean_depth_err*1e6:.2f} ppm")
    print(f"CO2-band contrast (this script) = {co2_excess_ppm:.1f} ppm ({co2_sigma:.1f} sigma band S/N)")
    print(f"H2O/H2S-band contrast (this script) = {h2o_excess_ppm:.1f} ppm ({h2o_sigma:.1f} sigma band S/N)")
    print(f"Fu et al. (2024) retrieval significances: H2O {PAPER_SIGMA['H2O']} sigma, CO2 {PAPER_SIGMA['CO2']} sigma, CO {PAPER_SIGMA['CO']} sigma, H2S {PAPER_SIGMA['H2S']} sigma")
    for family in ("free", "grid"):
        full = next(row for row in model_rows if row["family"] == family and str(row["model"]).startswith("FULL"))
        print(f"{family} full median: chi2/dof={float(full['chi2']):.2f}/{int(full['dof'])}")


if __name__ == "__main__":
    main()
