#!/usr/bin/env python3
"""Validate distributed files and reconstruct cluster endpoints (standard library only)."""
import csv
import hashlib
import math
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def read(name):
    with (ROOT / name).open(newline="") as handle:
        return list(csv.DictReader(handle))

# Every release file has a recorded digest; reruns may intentionally change outputs.
for row in read("MANIFEST.csv"):
    path = ROOT / row["path"]
    assert path.is_file(), row["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], row["path"]
plots = read("data/clean_plot_data_public.csv")
organs = read("data/cluster_organ_data_public.csv")
programs = read("data/retained_program_codebook_public.csv")
assert len(plots) == 264 and len(programs) == 66
assert len({r["Public_Plot_ID"] for r in plots}) == 264
ids = {r["Public_Program_ID"] for r in programs}
assert len(ids) == 66 and {r["Public_Program_ID"] for r in plots} == ids
by_plot = defaultdict(list)
for row in organs:
    assert row["Public_Program_ID"] in ids
    by_plot[row["Public_Plot_ID"]].append(row)
assert set(by_plot) == {r["Public_Plot_ID"] for r in plots}
for row in plots:
    observations = by_plot[row["Public_Plot_ID"]]
    assert len(observations) == 20
    assert all(x["Public_Program_ID"] == row["Public_Program_ID"] for x in observations)
    values = [float(x["Analysis_Rating_Pct"]) for x in observations if x["Analysis_Rating_Pct"] not in ("", "NA")]
    assert all(0 <= x <= 100 for x in values)
    assert len(values) == int(row["N_Clusters_Rated"])
    positive = sum(x > 0 for x in values)
    assert positive == int(row["Cluster_Positive_Count"])
    assert math.isclose(100 * positive / len(values), float(row["Cluster_Incidence_Pct"]), abs_tol=1e-10)
    assert math.isclose(sum(values) / len(values), float(row["Cluster_Severity_Pct"]), abs_tol=1e-10)
# Check Figure 6 overlap rules independently of the plotting implementation.
events = read("descriptive_newa_timeline/NEWA_2021_2025_black_rot_grouped_periods.csv")
counts = []
for row in read("data/annual_methods_weather.csv"):
    start = datetime.fromisoformat(row["Experimental_Application_Dates"].split("; ")[0])
    end = datetime.fromisoformat(row["Final_Assessment_Date"]) + timedelta(days=1)
    counts.append(sum(e["year"] == row["Year"] and datetime.fromisoformat(e["end"]) >= start and datetime.fromisoformat(e["start"]) < end for e in events))
assert counts == [21, 10, 10, 14, 13], counts
print("PASS: file checksums; 264 plots; 66 program-years; all cluster endpoints; timeline totals", counts)
