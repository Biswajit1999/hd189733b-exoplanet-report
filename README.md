# HD 189733 b — Exoplanet Atmosphere Report

The "blue planet" — famous for scattered-light Rayleigh haze that gives it
its inferred azure color — and one of the most intensively studied hot
Jupiters ever found. This repo re-derives a real, high-precision JWST NIRCam
transmission spectrum with molecular-band detections at extremely high
statistical significance.

**[Open the full report](index.html)** (open locally in a browser, or serve
with `python -m http.server` from this directory).

## What's real here

- **System parameters** — queried live from the NASA Exoplanet Archive TAP
  service (`pscomppars` table).
- **JWST NIRCam spectrum** — 139 real wavelength bins, 2.4-5.0 microns,
  from a study reporting H2S and metal enrichment in this atmosphere,
  released publicly on Zenodo
  ([10.5281/zenodo.11459715](https://doi.org/10.5281/zenodo.11459715)).
- **Analysis** — `scripts/analyze_spectrum.py` computes the weighted mean
  transit depth and compares two real molecular-absorption windows (CO2 near
  4.3-4.6 um, H2O/H2S-influenced near 2.6-3.0 um) against a nearby
  featureless continuum, with a real statistical significance for each. Run
  it yourself:

  ```bash
  pip install -r requirements.txt
  python scripts/analyze_spectrum.py
  ```

## Repository structure

```text
index.html              the report webpage
data/                    real JWST NIRCam transmission spectrum (Zenodo)
scripts/analyze_spectrum.py   real band-comparison analysis
figures/                 generated plot + summary_statistics.csv
```

## Key finding this repo shows directly

Weighted mean transit depth 24183 ppm. Real CO2-band excess: 239 ppm at
22.1-sigma significance. Real H2O/H2S-band excess: 255 ppm at 24.7-sigma
significance. Both are firm, high-confidence detections in this specific
data — a useful contrast against this series' K2-18 b report, where an
identically structured band-vs-continuum comparison came out consistent with
zero (~0.2 sigma), directly illustrating the real difference between a
well-characterized bright target and a marginal, disputed one.

## References

1. Bouchy, F. et al., 2005. A very hot Jupiter transiting the bright K star
   HD 189733. *Astronomy & Astrophysics*, 444(1), pp.L15-L19.
2. Knutson, H.A. et al., 2007. A map of the day-night contrast of the
   extrasolar planet HD 189733b. *Nature*, 447, pp.183-186.
3. Pont, F. et al., 2013. The prevalence of dust on the exoplanet HD
   189733b from Hubble Space Telescope camera spectroscopy. *Monthly
   Notices of the Royal Astronomical Society*, 432(4), pp.2917-2944.
4. Zenodo record
   [10.5281/zenodo.11459715](https://doi.org/10.5281/zenodo.11459715),
   "Products for 'Hydrogen sulfide and metal-enriched atmosphere for a
   Jupiter-mass exoplanet.'"
5. NASA Exoplanet Archive, <https://exoplanetarchive.ipac.caltech.edu/>.

## Author

Biswajit Jana — [Portfolio](https://biswajit1999.github.io/Biswajit_Jana.github.io/) · [GitHub](https://github.com/Biswajit1999) · [LinkedIn](https://www.linkedin.com/in/biswajit-jana-27011a151/) · [ORCID](https://orcid.org/0009-0002-2411-1891)
