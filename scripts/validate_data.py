"""Validate the Zenodo v1 spectrum/model subset committed here."""

import hashlib
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
EXPECTED_SET_SHA256 = "8e387ba124007b43e3123e91ff5109c0ed0e4aa2c5a1162e861f666c0b284876"


def main() -> None:
    paths = [DATA_DIR / "nircam_transmission_spectrum.txt"] + sorted(
        (DATA_DIR / "models").rglob("*.txt")
    )
    digest = hashlib.sha256()
    for path in paths:
        relative = path.relative_to(DATA_DIR).as_posix()
        text = path.read_text(encoding="utf-8").replace("\r\n", "\n").replace("\r", "\n")
        digest.update(relative.encode("utf-8") + b"\0" + text.encode("utf-8"))
        table = np.loadtxt(path)
        if table.ndim != 2 or table.shape[1] not in {2, 4}:
            raise SystemExit(f"unexpected table shape for {relative}: {table.shape}")
    if digest.hexdigest() != EXPECTED_SET_SHA256:
        raise SystemExit(f"data-set SHA-256 mismatch: {digest.hexdigest()}")
    print(f"Validated spectrum plus {len(paths) - 1} supplied model products")


if __name__ == "__main__":
    main()
