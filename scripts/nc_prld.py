"""Nucleocapsid knuckle against prion-like domain, on the HMM knuckle calls.

Re-runs every CCHC-based result on `nc_hits.tsv` from scripts/nc_scan.py. The
earlier versions of these numbers used a regex and are superseded -- it called a
knuckle in 70% of copia Gags where the profile calls one in 21%.

Three definitions are reported side by side, because the answer moves with the
definition and that fact is part of the result:

    strict     Pfam gathering threshold -- 61% sensitive against InterPro
    sensitive  Pfam E < 1e-3            -- 87%
    combined   Pfam E < 1e-3 or PROSITE -- 92% at 99% specificity

`combined` is the primary. It is the closest local approximation to InterProScan,
calibrated on the 336 GyDB cores where InterPro calls are known. A definition
that misses knuckles does not bias the odds ratio randomly -- it attenuates it
toward 1, which is how an earlier pass on `strict` alone concluded there was no
effect when there was one.

Length is in every model because it dominates everything here: a prion-like
domain is roughly 15x more likely per e-fold of Gag length, and any feature that
correlates with length will look associated until length is controlled.

Caveat that no model here addresses: these elements are related. The odds ratios
treat thousands of tips of one phylogeny as independent draws.

Needs the conda environment (scipy). Run from the repo root:
  python scripts/nc_prld.py
"""

import csv
import math
import os
import statistics as st
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

FULL = os.path.join("data", "gag_full")
CLASSIFICATION = os.path.join("data", "REXdb", "Viridiplantae_v4.0.classification")


def read_tsv(path):
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def read_fasta(path):
    header, parts = None, []
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if line.startswith(">"):
            if header:
                yield header, "".join(parts)
            header, parts = line[1:], []
        elif line:
            parts.append(line)
    if header:
        yield header, "".join(parts)


def logistic(rows, key):
    """PrLD ~ knuckle + log(length). Returns (OR, lo, hi, p) or None."""
    y = np.array([float(r["prd"]) for r in rows])
    x = np.array([float(r[key]) for r in rows])
    if y.sum() < 5 or len(set(x)) < 2:
        return None
    X = np.column_stack([np.ones(len(y)), x,
                         np.array([math.log(r["len"]) for r in rows])])
    beta = np.zeros(3)
    for _ in range(100):
        p = 1 / (1 + np.exp(-(X @ beta)))
        W = np.clip(p * (1 - p), 1e-9, None)
        try:
            step = np.linalg.solve((X * W[:, None]).T @ X,
                                   (X * W[:, None]).T @ (X @ beta + (y - p) / W))
        except np.linalg.LinAlgError:
            return None
        if np.max(np.abs(step - beta)) < 1e-9:
            beta = step
            break
        beta = step
    p = 1 / (1 + np.exp(-(X @ beta)))
    try:
        cov = np.linalg.inv((X * np.clip(p * (1 - p), 1e-9, None)[:, None]).T @ X)
    except np.linalg.LinAlgError:
        return None
    se = np.sqrt(np.diag(cov))
    if abs(beta[1]) > 15 or se[1] > 15:
        return None
    return (math.exp(beta[1]), math.exp(beta[1] - 1.96 * se[1]),
            math.exp(beta[1] + 1.96 * se[1]),
            2 * (1 - stats.norm.cdf(abs(beta[1] / se[1]))))


def show(rows, label):
    cells = []
    for key in ("strict", "sensitive", "combined"):
        fit = logistic(rows, key)
        cells.append("OR %6.3f (%.3f-%.3f) p=%-9.2g" % fit if fit
                     else "%-34s" % "not estimable")
    print("   %-30s n=%5d PrLD=%4d  %s  %s  %s"
          % (label, len(rows), sum(r["prd"] for r in rows),
             cells[0], cells[1], cells[2]))


def main():
    lineage = {}
    for row in open(CLASSIFICATION, encoding="utf-8"):
        f = row.rstrip("\n").split("\t")
        if len(f) >= 5 and f[2] == "LTR":
            lineage[f[0]] = f[-1]

    nc = {r["rexdb_id"]: r for r in read_tsv(os.path.join(FULL, "nc_hits.tsv"))}
    plaac = {r["rexdb_id"]: r for r in read_tsv(os.path.join(FULL,
                                                             "gag_full_plaac.tsv"))}
    bounds = {r["rexdb_id"]: r for r in read_tsv(os.path.join(FULL, "gag_full.tsv"))}
    seqs = {}
    for short in ("copia", "gypsy"):
        seqs.update(dict(read_fasta(os.path.join(FULL, "gag_full_%s.faa" % short))))

    rows = []
    for element, seq in seqs.items():
        n, p, b = nc[element], plaac[element], bounds[element]
        rows.append({
            "id": element, "sf": b["superfamily"], "lin": lineage.get(element, "NA"),
            "len": len(seq), "clip": b["down_clipped_by"],
            "up": int(b["gag_upstream"]), "core": int(b["gag_core_len"]),
            "prd": p["has_prd"] == "1",
            "prd_start": int(p["prd_start"]) - 1 if p["prd_start"] else None,
            "prd_end": int(p["prd_end"]) if p["prd_end"] else None,
            "strict": n["nc_strict"] == "1", "sensitive": n["nc_sensitive"] == "1",
            "combined": n.get("nc_combined") == "1",
            "start": int(n["best_start"]) - 1 if n["best_start"] else None,
            "end": int(n["best_end"]) if n["best_end"] else None,
            "model": n["best_model"]})

    print("=== knuckle prevalence ===\n")
    print("%-11s %6s %11s %11s %11s"
          % ("", "n", "strict", "sensitive", "combined"))
    for sf in ("Ty1/copia", "Ty3/gypsy"):
        sub = [r for r in rows if r["sf"] == sf]
        cells = ["%5d %3.0f%%" % (sum(r[k] for r in sub),
                                  100 * sum(r[k] for r in sub) / len(sub))
                 for k in ("strict", "sensitive", "combined")]
        print("%-11s %6d %11s %11s %11s" % (sf, len(sub), *cells))
    print("\n   which model fired: %s"
          % ", ".join("%s %d" % kv for kv in
                      Counter(r["model"] for r in rows if r["strict"]).most_common()))

    print("\n=== by lineage (combined definition) ===\n")
    print("%-12s %6s %9s %9s %9s" % ("lineage", "n", "knuckle", "PrLD", "both"))
    for lin, n in Counter(r["lin"] for r in rows).most_common(12):
        sub = [r for r in rows if r["lin"] == lin]
        print("%-12s %6d %8.0f%% %8.1f%% %8.1f%%"
              % (lin, n, 100 * sum(r["combined"] for r in sub) / n,
                 100 * sum(r["prd"] for r in sub) / n,
                 100 * sum(r["combined"] and r["prd"] for r in sub) / n))

    print("\n=== where does the knuckle sit? (strict calls) ===\n")
    for sf in ("Ty1/copia", "Ty3/gypsy"):
        sub = [r for r in rows if r["sf"] == sf and r["strict"] and r["start"]]
        zone = Counter()
        for r in sub:
            zone["N-terminus" if r["start"] < r["up"] else
                 "capsid core" if r["start"] < r["up"] + r["core"] else
                 "C-terminal region"] += 1
        past = [r["start"] - (r["up"] + r["core"]) for r in sub]
        print("%-11s %5d knuckles: %s" % (sf, len(sub), dict(zone)))
        print("            median %+d aa past the capsid core (IQR %+d to %+d)"
              % (st.median(past), st.quantiles(past, n=4)[0],
                 st.quantiles(past, n=4)[2]))

    print("\n=== knuckle relative to the called domain ===\n")
    both = [r for r in rows if r["strict"] and r["prd"] and r["start"] is not None]
    rel = Counter()
    gaps = []
    for r in both:
        if r["prd_start"] <= r["start"] <= r["prd_end"]:
            rel["inside the domain"] += 1
        elif r["start"] > r["prd_end"]:
            rel["downstream of the domain"] += 1
            gaps.append(r["start"] - r["prd_end"])
        else:
            rel["upstream of the domain"] += 1
    print("   %d Gags carry both: %s" % (len(both), dict(rel)))
    if gaps:
        print("   median gap when downstream: %d aa" % st.median(gaps))

    print("\n=== the trade-off, both definitions ===")
    print("%-37s %-24s %-24s\n" % ("", "STRICT (Pfam GA)", "SENSITIVE (E<1e-3)"))
    cop = [r for r in rows if r["sf"] == "Ty1/copia"]
    gyp = [r for r in rows if r["sf"] == "Ty3/gypsy"]
    show(cop, "copia, all")
    show([r for r in cop if r["clip"] == "prot"], "copia, PROT-bounded only")
    show([r for r in cop if r["lin"] == "Ale"], "copia, Ale only")
    show([r for r in cop if r["lin"] == "Ivana"], "copia, Ivana only")
    show(gyp, "gypsy, all")
    show([r for r in gyp if r["lin"] not in
          {"Ogre", "TatI", "TatII", "TatIII", "Tatius", "Retand"}],
         "gypsy, excluding Tat + Retand")
    show([r for r in gyp if r["lin"] == "Athila"], "gypsy, Athila only")


if __name__ == "__main__":
    main()
