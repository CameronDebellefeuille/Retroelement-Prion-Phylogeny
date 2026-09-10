"""PLAAC over the entire GyDB cores database, grouped by protein domain.

The question: is Gag more prion-like than other retroelement proteins, or does
PLAAC hit everything at some background rate?

The confound that has to be handled is length. PLAAC's core window is 60 aa
(L-5), so a domain shorter than that cannot be called at all, and longer
proteins give the window more opportunity. Raw per-domain rates are therefore
compared only within length strata.

Run from the repo root:

  python scripts/plaac_all_cores.py                            # alpha 0.5, the ruler
  python scripts/plaac_all_cores.py <cores> <outdir> <alpha>

Writes the raw PLAAC table (all_cores_a05.tsv at alpha 0.5) next to its own
per-domain summary. Everything else PLAAC-derived is measured against that
raw table, so it must not be regenerated at a different alpha by accident.
"""

import math
import os
import statistics
import subprocess
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
PLAAC_JAR = os.path.join(REPO, "plaac", "plaac.jar")
GYDB = os.path.join(REPO, "data", "processed", "gydb")
RAW  = os.path.join(REPO, "data", "raw", "gydb")
CORE_LENGTH = 60
DEFAULT_ALPHA = "0.5"


def read_fasta(path):
    header, parts = None, []
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.strip()
        if line.startswith(">"):
            if header:
                yield header, "".join(parts)
            header, parts = line[1:].strip(), []
        elif line:
            parts.append(line)
    if header:
        yield header, "".join(parts)


def domain_of(header):
    """GAG_Ty1B -> GAG. Headers starting with '_' keep the leading part."""
    h = header.lstrip("_")
    return h.split("_")[0] if "_" in h else h


def run_plaac(fasta, alpha, raw_out=None):
    """Run the jar and, if asked, keep its raw table before parsing anything.

    Keeping the raw output is the whole point of raw_out. all_cores_a05.tsv --
    the alpha=0.5 run every LLR in the project is measured against -- used to
    have no producer at all, because this function captured stdout, pulled six
    columns out of it and dropped the other thirty-two on the floor.
    """
    result = subprocess.run(
        ["java", "-Xmx4g", "-jar", PLAAC_JAR, "-i", fasta,
         "-a", alpha, "-c", str(CORE_LENGTH)],
        capture_output=True, text=True, check=True)
    lines = [l for l in result.stdout.splitlines() if not l.startswith("#")]
    if raw_out:
        with open(raw_out, "w", encoding="utf-8", newline="\n") as handle:
            handle.write("\n".join(lines) + "\n")
        print("wrote %s (raw PLAAC output)" % raw_out)
    header = lines[0].split("\t")
    at = {c: header.index(c) for c in
          ("SEQid", "LLR", "NLLR", "PRDstart", "PRDend", "COREscore")}
    out = {}
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) < len(header):
            continue
        called = f[at["COREscore"]] not in ("NaN", "")
        out[f[at["SEQid"]]] = {
            "llr": f[at["LLR"]], "nllr": f[at["NLLR"]],
            "prd_start": f[at["PRDstart"]] if called else "",
            "prd_end": f[at["PRDend"]] if called else "",
            "has_prd": 1 if called else 0}
    return out


def two_proportion_p(k1, n1, k2, n2):
    """Two-sided z-test on two proportions. Plain normal approximation."""
    if n1 == 0 or n2 == 0:
        return float("nan")
    p1, p2 = k1 / n1, k2 / n2
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    if se == 0:
        return float("nan")
    z = (p1 - p2) / se
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


def main():
    cores = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RAW, "cores-database")
    outdir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(GYDB, "plaac_data")
    alpha = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_ALPHA
    os.makedirs(outdir, exist_ok=True)
    raw_out = os.path.join(outdir, "all_cores_a05.tsv" if alpha == "0.5"
                           else "all_cores_a%s.tsv" % alpha.replace(".", ""))

    seqs = dict(read_fasta(cores))
    print("%d sequences in %s" % (len(seqs), cores))

    fasta = os.path.join(outdir, "all_cores.faa")
    with open(fasta, "w", encoding="utf-8", newline="\n") as handle:
        for h, s in seqs.items():
            handle.write(">%s\n%s\n" % (h, s))

    print("running plaac (alpha=%s, core=%d) ..." % (alpha, CORE_LENGTH), flush=True)
    plaac = run_plaac(fasta, alpha, raw_out=raw_out)

    rows = []
    for h, s in seqs.items():
        r = plaac.get(h, {"llr": "", "nllr": "", "prd_start": "",
                          "prd_end": "", "has_prd": 0})
        rows.append({"header": h, "domain": domain_of(h), "length": len(s), **r})

    path = os.path.join(outdir, "all_cores_plaac.tsv")
    cols = ["header", "domain", "length", "llr", "nllr",
            "prd_start", "prd_end", "has_prd"]
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["domain"], r["header"])):
            handle.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("wrote %s\n" % path)

    by = defaultdict(list)
    for r in rows:
        by[r["domain"]].append(r)

    print("=== raw rate by domain (n >= 15), NOT length-controlled ===")
    print("%-12s %6s %8s %7s %8s" % ("domain", "n", "medlen", "PRD", "rate"))
    for d in sorted(by, key=lambda d: -sum(r["has_prd"] for r in by[d]) / len(by[d])):
        g = by[d]
        if len(g) < 15:
            continue
        k = sum(r["has_prd"] for r in g)
        print("%-12s %6d %8d %7d %7.1f%%"
              % (d, len(g), statistics.median(r["length"] for r in g), k,
                 100 * k / len(g)))

    scorable = [r for r in rows if r["length"] >= CORE_LENGTH]
    print("\nsequences shorter than the %d aa core window (cannot be called): %d"
          % (CORE_LENGTH, len(rows) - len(scorable)))

    print("\n=== GAG vs everything else, within length strata ===")
    print("%-14s %18s %18s %10s" % ("length band", "GAG", "other", "p"))
    bands = [(60, 200), (200, 300), (300, 400), (400, 600), (600, 10000)]
    tg = to = kg = ko = 0
    for lo, hi in bands:
        g = [r for r in scorable if r["domain"] == "GAG" and lo <= r["length"] < hi]
        o = [r for r in scorable if r["domain"] != "GAG" and lo <= r["length"] < hi]
        if not g or not o:
            continue
        k1, k2 = sum(r["has_prd"] for r in g), sum(r["has_prd"] for r in o)
        tg += len(g); to += len(o); kg += k1; ko += k2
        print("%-14s %7d/%-4d %5.1f%% %7d/%-4d %5.1f%% %10.3g"
              % ("%d-%d" % (lo, hi), k1, len(g), 100 * k1 / len(g),
                 k2, len(o), 100 * k2 / len(o), two_proportion_p(k1, len(g), k2, len(o))))
    print("%-14s %7d/%-4d %5.1f%% %7d/%-4d %5.1f%% %10.3g"
          % ("ALL >=60", kg, tg, 100 * kg / tg, ko, to, 100 * ko / to,
             two_proportion_p(kg, tg, ko, to)))


if __name__ == "__main__":
    main()
