# Data source

`nircam_transmission_spectrum.txt` is downloaded, unmodified, from Zenodo
record **10.5281/zenodo.11459715** ("Products for 'Hydrogen sulfide and
metal-enriched atmosphere for a Jupiter-mass exoplanet'"), file
`Data/transit spectrum/HD189733b_NIRCam_spectrum.txt`.

Retrieved: 2026-08-11, via `https://zenodo.org/api/records/11459715`.

Four whitespace-separated columns, one header line (`#`-prefixed):

1. wavelength [micron]
2. wavelength bin half-width [micron]
3. transit depth, (Rp/Rs)^2
4. 1-sigma uncertainty on the transit depth

The repository also preserves all 12 median model curves under `models/free/`
and `models/grid/` from the same v1 archive. `FULL` is the supplied full model;
`xMOL` curves remove the named absorber, and `Cloud_only` retains only the
cloud contribution. These are publication products, not fits performed here.

The downloaded archive MD5 is `7b39f9b6fdfc0d19a919a6d27eebde28`.
The LF-normalized, filename-delimited SHA-256 over the committed spectrum and
12 model files is `8e387ba124007b43e3123e91ff5109c0ed0e4aa2c5a1162e861f666c0b284876`.
Run `python scripts/validate_data.py` to verify the subset.
