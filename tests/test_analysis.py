"""Executable checks on the weighted-mean/band-contrast statistics and
a regression guard that the pipeline still reproduces the documented
headline numbers when run on the real downloaded data."""

import csv

import numpy as np
import analyze_spectrum as spec


def test_weighted_mean_matches_hand_computed_case():
    values = np.array([1.0, 2.0])
    errors = np.array([1.0, 0.5])  # weights 1 and 4
    mean, err = spec.weighted_mean(values, errors)
    assert np.isclose(mean, 1.8, rtol=1e-10)
    assert np.isclose(err, np.sqrt(1.0 / 5.0), rtol=1e-10)


def test_band_stats_selects_only_points_in_range():
    wave = np.array([1.0, 2.0, 3.0, 4.0])
    depth = np.array([10.0, 20.0, 30.0, 40.0])
    err = np.full(4, 1.0)
    mean, _ = spec.band_stats(wave, depth, err, (1.5, 3.5))
    # Only the two points at wave=2.0 and 3.0 fall in range -> plain mean 25.
    assert np.isclose(mean, 25.0)


def test_pipeline_reproduces_documented_headline_numbers():
    spec.FIG_DIR.mkdir(exist_ok=True)
    spec.main()
    rows = {}
    with (spec.FIG_DIR / "summary_statistics.csv").open() as f:
        for row in csv.DictReader(f):
            rows[row["quantity"]] = row["value"]
    assert int(rows["n_wavelength_bins"]) == 139
    assert abs(float(rows["co2_band_contrast_snr_this_script"]) - 22.1) < 0.1
    assert abs(float(rows["h2o_h2s_band_contrast_snr_this_script"]) - 24.7) < 0.1
    # The paper's own real retrieval significances must stay attached for
    # comparison, not silently dropped.
    assert float(rows["paper_retrieval_significance_H2O"]) == 13.4
    assert float(rows["paper_retrieval_significance_CO2"]) == 11.2
