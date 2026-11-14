"""Phase 1 audit helpers. No writes outside the new revision directory."""
from pathlib import Path
import csv
import hashlib
import json
import math

REV = Path(__file__).resolve().parents[1]
PROJECT = REV.parent
P = "XenoPairDB_Integrity_Pipeline"
SUB = "1.投稿/5.Xenotransplanation"
SNAP = REV / "inputs/snapshot"
NA = {"", "na", "n/a", "nan", "none", "null"}

def number(value):
    if value is None or str(value).strip().lower() in NA:
        return None
    x = float(value)
    if not math.isfinite(x):
        raise ValueError(f"Non-finite number: {value}")
    return x

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else [])
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: "NA" if v is None else v for k, v in r.items()})

def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")

def canonical(name):
    return SNAP / P / "canonical_links" / name

def unique(rows, field="pair_id"):
    seen = set()
    for r in rows:
        if not r[field] or r[field] in seen:
            raise ValueError(f"Missing or duplicate {field}: {r[field]}")
        seen.add(r[field])

def identity_key(r):
    # Deliberately conservative: no blind alias collapsing or endpoint merging.
    return json.dumps([r[k].strip() for k in
                       ["pig_gene", "human_partner", "direction", "phenotype"]],
                      ensure_ascii=False, separators=(",", ":"))

def stable_id(r):
    return "XPA-" + hashlib.sha256(identity_key(r).encode()).hexdigest()[:20]

def short_donor(name):
    return name.split(" (")[0]
