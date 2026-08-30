"""XenoPairDB: a quantitative resource of molecular incompatibilities in
pig-to-human xenotransplantation.

Submodules
----------
scoring       : ISI v2, Tier-2 atlas score, antigenicity v2 (pure functions)
benchmark     : seven-model (M0-M5) resampling + LOSO benchmark framework
antigenicity  : surfaceome antigenicity screen and external literature benchmark
"""

from . import antigenicity, benchmark, scoring

__version__ = "1.0.0"
__all__ = ["scoring", "benchmark", "antigenicity", "data_path"]

from pathlib import Path as _Path


def data_path(data_dir, filename):
    """Resolve ``filename`` inside ``data_dir``, searching sub-directories.

    The Zenodo deposit ships the tables in a layered layout
    (``tier1_gold_standard/``, ``tier2_atlas/``, ``tier3_surfaceome/``,
    ``benchmark/``, ``clinical/``).  Pointing ``--data-dir`` at the unpacked
    deposit root therefore has to work without the user first flattening the
    files by hand; a flat directory keeps working because the direct hit is
    tried first.
    """
    root = _Path(data_dir)
    direct = root / filename
    if direct.is_file():
        return str(direct)
    matches = sorted(q for q in root.rglob(filename) if q.is_file())
    if len(matches) == 1:
        return str(matches[0])
    if len(matches) > 1:
        raise FileNotFoundError(
            f"{filename!r} is ambiguous under {root}: "
            + ", ".join(str(m) for m in matches))
    subdirs = sorted({str(q) for q in root.rglob("*") if q.is_dir()})
    raise FileNotFoundError(
        f"{filename!r} not found under {root}. Searched: "
        + (", ".join([str(root)] + subdirs) if subdirs else str(root)))

