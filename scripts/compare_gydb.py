"""Re-test the GyDB findings on REXdb, section by section.

GYDB_FINDINGS.md makes five claims that REXdb can check independently. GyDB is
2,636 protein cores across four superfamilies, cut per element with no common
rule; REXdb here is 13,572 plant Gags cut by one rule. If a claim survives both,
it is unlikely to be an artefact of either database's boundaries.

Each section prints the GyDB number beside the REXdb one, so a claim that fails
fails visibly rather than quietly.

Everything is reported twice where it matters: over all Gags, and over the 1,121
cluster representatives that are the tree tips. The representatives are not
independent either -- they are related tips of one phylogeny -- but 80% identity
clustering removes the grossest duplication, and a claim that only holds on the
redundant set is a claim about which elements got sequenced.

Needs the conda environment (scipy). Run from the repo root:
  python scripts/compare_gydb.py
"""

import csv
import math
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

FULL = os.path.join("data", "gag_full")
PLAAC = os.path.join(FULL, "gag_full_plaac.tsv")
TRAITS = os.path.join(FULL, "traits_full.tsv")
DOMAINS = os.path.join("data", "domain_scan", "rexdb_domains_plaac.tsv")
CLUSTER = os.path.join("data", "cluster")
REPORT = os.path.join(FULL, "gydb_check.txt")

CCHC = re.compile(r"C.{2}C.{4}H.{4}C")
BANDS = [(0, 200), (200, 400), (400, 600), (600, 10 ** 6)]


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


def read_tsv(path):
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def logistic(y, X, names):
    """Plain IRLS, so the odds ratios are comparable to GyDB's logistic step."""
    X = np.column_stack([np.ones(len(y)), X])
    beta = np.zeros(X.shape[1])
    for _ in range(50):
        eta = X @ beta
        p = 1 / (1 + np.exp(-eta))
        W = np.clip(p * (1 - p), 1e-9, None)
        z = eta + (y - p) / W
        beta_new = np.linalg.solve((X * W[:, None]).T @ X, (X * W[:, None]).T @ z)
        if np.max(np.abs(beta_new - beta)) < 1e-8:
            beta = beta_new
            break
        beta = beta_new
    p = 1 / (1 + np.exp(-(X @ beta)))
    cov = np.linalg.inv((X * np.clip(p * (1 - p), 1e-9, None)[:, None]).T @ X)
    se = np.sqrt(np.diag(cov))
    out = []
    for i, name in enumerate(["intercept"] + names):
        if name == "intercept":
            continue
        lo, hi = beta[i] - 1.96 * se[i], beta[i] + 1.96 * se[i]
        pval = 2 * (1 - stats.norm.cdf(abs(beta[i] / se[i])))
        out.append((name, math.exp(beta[i]), math.exp(lo), math.exp(hi), pval))
    return out


def mantel_haenszel(strata):
    """Common odds ratio across strata of (exposed_case, exposed_ctl, unexp_case, unexp_ctl)."""
    num = den = 0.0
    for a, b, c, d in strata:
        n = a + b + c + d
        if n == 0:
            continue
        num += a * d / n
        den += b * c / n
    return num / den if den else float("nan")


def section(title, gydb):
    print("\n" + "=" * 78)
    print(title)
    print("GyDB: %s" % gydb)
    print("-" * 78)


def main():
    plaac = {r["rexdb_id"]: r for r in read_tsv(PLAAC)}
    traits = {r["rexdb_id"]: r for r in read_tsv(TRAITS)}
    seqs = {}
    for short in ("copia", "gypsy"):
        seqs.update(dict(read_fasta(os.path.join(FULL,
                                                 "gag_full_%s.faa" % short))))
    tips = set()
    for short in ("copia", "gypsy"):
        tips.update(r["rexdb_id"] for r in
                    read_tsv(os.path.join(CLUSTER, "tips_%s.tsv" % short)))
    print("REXdb: %d full-length Gags, %d of them cluster representatives"
          % (len(seqs), len(tips & set(seqs))))

    # ---- 3.1 is prion-like character confined to Gag? ----------------------
    section("3.1  Prion-like character is confined to Gag",
            "26/26 called domains in GAG; 0 in 2,172 non-Gag proteins")
    by_domain = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for row in read_tsv(DOMAINS):
        length = int(row["length"])
        band = next(b for b in BANDS if b[0] <= length < b[1])
        cell = by_domain[row["domain"]][band]
        cell[0] += int(row["has_prd"])
        cell[1] += 1
    print("%-10s %s" % ("domain", "  ".join("%12s" % ("%d-%d" % b if b[1] < 10 ** 6
                                                      else "%d+" % b[0])
                                            for b in BANDS)))
    for domain in sorted(by_domain, key=lambda d: -sum(v[1] for v in
                                                       by_domain[d].values())):
        cells = []
        for band in BANDS:
            hit, n = by_domain[domain][band]
            cells.append("%12s" % ("%d/%d" % (hit, n) if n else "-"))
        print("%-10s %s" % (domain, "  ".join(cells)))

    # ---- 3.2 where in Gag does the domain sit? -----------------------------
    section("3.2  Domain position splits by superfamily",
            "gypsy 0 N-terminal of 17; copia 2 N-terminal of 3")
    # The two databases do not label position the same way. REXdb's prd_region
    # is a two-way split at the capsid core boundary -- the boundary this
    # project cares about -- while GyDB used terciles of the protein
    # (<33% N-terminal, 33-67% middle, >67% C-terminal). Both are reported, or
    # the comparison is between two different questions: 53 of the 285 REXdb
    # calls that are N-terminal-of-the-core sit past 33% of the protein.
    print("REXdb rule: N-term = domain midpoint upstream of the capsid core\n")
    for name, keep in (("all Gags", set(seqs)), ("tree tips", tips & set(seqs))):
        counts = defaultdict(Counter)
        totals = Counter()
        for element in keep:
            row = plaac[element]
            totals[row["superfamily"]] += 1
            if row["has_prd"] == "1":
                counts[row["superfamily"]][row["prd_region"]] += 1
        for superfamily in sorted(totals):
            c = counts[superfamily]
            total = sum(c.values())
            print("%-9s %-10s n=%5d  PrLD %3d (%4.1f%%)  N-term %3d  middle %3d  C-term %3d"
                  % (name, superfamily, totals[superfamily], total,
                     100 * total / totals[superfamily],
                     c["N-term"], c["middle"], c["C-term"]))

    print("\nGyDB rule: terciles of the protein, so the two are comparable\n")
    for name, keep in (("all Gags", set(seqs)), ("tree tips", tips & set(seqs))):
        counts = defaultdict(Counter)
        totals = Counter()
        for element in keep:
            row = plaac[element]
            totals[row["superfamily"]] += 1
            if row["has_prd"] == "1" and row["prd_start"]:
                mid = (int(row["prd_start"]) + int(row["prd_end"])) / 2
                frac = mid / int(row["gag_len"])
                zone = ("N-term" if frac < 1 / 3 else
                        "middle" if frac <= 2 / 3 else "C-term")
                counts[row["superfamily"]][zone] += 1
        for superfamily in sorted(totals):
            c = counts[superfamily]
            total = sum(c.values())
            print("%-9s %-10s n=%5d  PrLD %3d (%4.1f%%)  N-term %3d  middle %3d  C-term %3d"
                  % (name, superfamily, totals[superfamily], total,
                     100 * total / totals[superfamily],
                     c["N-term"], c["middle"], c["C-term"]))

    # ---- 3.3 composition of the called domains -----------------------------
    section("3.3  PLAAC detects Q/N in the absence of charge",
            "PrLD Q+N 31.3% vs rest 9.1%; D+E 1.7% vs 11.2%; K+R 9.1% vs 12.9%")
    inside, outside = Counter(), Counter()
    for element, seq in seqs.items():
        row = plaac[element]
        if row["has_prd"] == "1" and row["prd_start"]:
            start, end = int(row["prd_start"]) - 1, int(row["prd_end"])
            inside.update(seq[start:end])
            outside.update(seq[:start] + seq[end:])
        else:
            outside.update(seq)
    n_in, n_out = sum(inside.values()), sum(outside.values())
    print("%d residues inside called domains, %d outside" % (n_in, n_out))
    print("%-6s %8s %8s %8s" % ("", "PrLD", "rest", "diff"))
    for label, group in (("Q", "Q"), ("N", "N"), ("P", "P"), ("G", "G"),
                         ("Y", "Y"), ("Q+N", "QN"), ("D+E", "DE"), ("K+R", "KR")):
        a = 100 * sum(inside[c] for c in group) / n_in
        b = 100 * sum(outside[c] for c in group) / n_out
        print("%-6s %7.1f%% %7.1f%% %+7.1f" % (label, a, b, a - b))

    # ---- 3.5 zinc knuckles vs prion-like domains ---------------------------
    section("3.5  Zinc knuckles and PrLDs are anti-correlated",
            "CCHC in 58% of Gags without a PrLD, 20% of those with; "
            "logistic OR 0.100 (0.032-0.312) adjusting for length")
    for name, keep in (("all Gags", sorted(seqs)), ("tree tips", sorted(tips & set(seqs)))):
        y, cchc, loglen, gypsy = [], [], [], []
        for element in keep:
            seq = seqs[element]
            y.append(int(plaac[element]["has_prd"] == "1"))
            cchc.append(int(bool(CCHC.search(seq))))
            loglen.append(math.log(len(seq)))
            gypsy.append(int(plaac[element]["superfamily"] == "Ty3/gypsy"))
        y = np.array(y); cchc = np.array(cchc)
        loglen = np.array(loglen); gypsy = np.array(gypsy)
        a = int(((y == 1) & (cchc == 1)).sum()); b = int(((y == 1) & (cchc == 0)).sum())
        c = int(((y == 0) & (cchc == 1)).sum()); d = int(((y == 0) & (cchc == 0)).sum())
        odds, fisher_p = stats.fisher_exact([[a, b], [c, d]])
        print("\n%s (n=%d): CCHC in %.0f%% of PrLD-bearing Gags, %.0f%% of the rest"
              % (name, len(y), 100 * a / max(a + b, 1), 100 * c / max(c + d, 1)))
        print("   raw 2x2 OR %.3f, Fisher p = %.2g" % (odds, fisher_p))

        # length quartiles, because length predicts PrLD strongly in both databases
        edges = np.quantile(np.exp(loglen), [0.25, 0.5, 0.75])
        strata = []
        for lo, hi in zip([0] + list(edges), list(edges) + [10 ** 9]):
            m = (np.exp(loglen) >= lo) & (np.exp(loglen) < hi)
            strata.append((int(((y == 1) & (cchc == 1) & m).sum()),
                           int(((y == 0) & (cchc == 1) & m).sum()),
                           int(((y == 1) & (cchc == 0) & m).sum()),
                           int(((y == 0) & (cchc == 0) & m).sum())))
        print("   Mantel-Haenszel common OR across length quartiles: %.3f"
              % mantel_haenszel(strata))
        for term, orr, lo, hi, p in logistic(
                y, np.column_stack([cchc, loglen, gypsy]),
                ["CCHC", "log(len)", "gypsy"]):
            print("   logistic %-10s OR %7.3f (%.3f-%.3f)  p = %.2g"
                  % (term, orr, lo, hi, p))
        med_prd = np.median(np.exp(loglen)[y == 1]) if (y == 1).any() else float("nan")
        print("   median Gag length: PrLD %.0f aa, no PrLD %.0f aa (Mann-Whitney p = %.2g)"
              % (med_prd, np.median(np.exp(loglen)[y == 0]),
                 stats.mannwhitneyu(np.exp(loglen)[y == 1],
                                    np.exp(loglen)[y == 0]).pvalue))

    # ---- 3.6 does disorder travel with prion-likeness? ---------------------
    section("3.6  Disorder and prion-likeness travel together",
            "called 0.50 of core disordered vs 0.29 for the rest")
    for name, keep in (("all Gags", set(seqs)), ("tree tips", tips & set(seqs))):
        hit, miss = [], []
        for element in keep:
            row = traits.get(element)
            if not row or row["disorder_scored"] != "1":
                continue
            (hit if plaac[element]["has_prd"] == "1" else miss).append(
                float(row["frac_disordered"]))
        print("%-9s PrLD %.2f (n=%d) vs no PrLD %.2f (n=%d), Mann-Whitney p = %.2g"
              % (name, np.median(hit), len(hit), np.median(miss), len(miss),
                 stats.mannwhitneyu(hit, miss).pvalue))


if __name__ == "__main__":
    with open(REPORT, "w", encoding="utf-8", newline="\n") as handle:
        class Tee:
            def write(self, text):
                sys.__stdout__.write(text)
                handle.write(text)

            def flush(self):
                sys.__stdout__.flush()
        sys.stdout = Tee()
        main()
        sys.stdout = sys.__stdout__
    print("\nwrote %s" % REPORT)
