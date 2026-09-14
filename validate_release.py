from pathlib import Path
import csv
import math
import sys

ROOT = Path(__file__).resolve().parents[1]
forbidden_names = {"subject_id", "hadm_id", "stay_id", "ecg_id", "charttime"}
errors = []

for p in ROOT.rglob("*.csv"):
    with p.open(newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader, [])
    found = forbidden_names.intersection(x.lower() for x in header)
    if found:
        errors.append(f"{p.relative_to(ROOT)} contains restricted identifier columns: {sorted(found)}")

required = [
    ROOT / "results" / "phase3b_hc3_precision_results.csv",
    ROOT / "results" / "phase3b_replication_check.csv",
    ROOT / "results" / "drug_flow_counts.csv",
    ROOT / "README.md",
    ROOT / "CITATION.cff",
    ROOT / "figshare_metadata.json",
]
for p in required:
    if not p.exists():
        errors.append(f"Missing required file: {p.relative_to(ROOT)}")

if not errors:
    with (ROOT / "results" / "phase3b_hc3_precision_results.csv").open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    by_drug = {r["ingredient"]: r for r in rows}
    checks = {
        "amiodarone": (6.93, 663),
        "propofol": (3.06, 2036),
        "haloperidol": (2.49, 662),
    }
    for drug, (estimate, n) in checks.items():
        r = by_drug.get(drug)
        if r is None:
            errors.append(f"Missing result: {drug}")
            continue
        if abs(float(r["hc3_adjusted_ms"]) - estimate) > 0.01:
            errors.append(f"Unexpected HC3 estimate for {drug}")
        if int(float(r["complete_cases"])) != n:
            errors.append(f"Unexpected complete-case count for {drug}")

if errors:
    print("RELEASE VALIDATION FAILED")
    for e in errors:
        print("-", e)
    sys.exit(1)

print("RELEASE VALIDATION PASSED")
print("No public CSV contains restricted identifier columns.")
