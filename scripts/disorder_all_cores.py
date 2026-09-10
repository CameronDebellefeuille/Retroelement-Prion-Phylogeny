"""metapredict over every GyDB core, so disorder can be compared by protein class.

The companion to the PLAAC run. Unlike PLAAC there is no 60 aa window, so the
short classes -- CHR especially, median 50 aa -- can be scored here even though
they had to be dropped from the PLAAC figure. The n per class therefore differs
between the two figures, deliberately.

PLAAC calls are read straight from its own output (all_cores_a05.tsv, alpha=0.5)
rather than from a reshaped intermediate. Only the PrLD call is taken from it;
protein class is derived here, as it always was.

Run under WSL in the copia-prld env:
  python scripts/disorder_all_cores.py
"""

import csv
import os

import metapredict as meta

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
GYDB = os.path.join(REPO, "data", "processed", "gydb")
RAW  = os.path.join(REPO, "data", "raw", "gydb")
CORES = os.path.join(RAW, "cores-database")
TWO = {"PR_ULP1", "PR_OTU", "ORF1_Nterdomain"}


def protein_class(header):
    parts = header.split("_")
    two = "_".join(parts[:2])
    return two if two in TWO else parts[0]


def read_cores(path):
    header, parts = None, []
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.strip()
        if line.startswith(">"):
            if header:
                yield header, "".join(parts)
            header, parts = line[1:], []
        elif line:
            parts.append(line)
    if header:
        yield header, "".join(parts)


def main():
    # PLAAC writes PRDlen 0 when nothing was called. One core name carries a
    # trailing space in its output, which has to go or the join on header fails.
    plaac = {r["SEQid"].strip(): ("1" if int(r["PRDlen"]) > 0 else "0")
             for r in csv.DictReader(open(os.path.join(GYDB, "plaac_data", "all_cores_a05.tsv"),
                                          encoding="utf-8"), delimiter="\t")
             if r.get("SEQid")}
    # metapredict refuses anything outside the 20 standard residues. Substituting
    # a residue for X would invent sequence, so those cores are dropped and
    # counted instead.
    standard = set("ACDEFGHIKLMNPQRSTVWY")
    rows, short, nonstandard = [], 0, 0
    for header, seq in read_cores(CORES):
        clean = "".join(c for c in seq.upper() if c.isalpha())
        if len(clean) < 20:            # metapredict needs something to work with
            short += 1
            continue
        if not set(clean) <= standard:
            nonstandard += 1
            continue
        result = meta.predict_disorder_domains(clean)
        in_idr = sum(e - s for s, e in result.disordered_domain_boundaries)
        rows.append({
            "seqid": header,
            "protein_class": protein_class(header),
            "prot_len": len(clean),
            "mean_disorder": round(sum(result.disorder) / len(result.disorder), 4),
            "frac_disordered": round(in_idr / len(clean), 4),
            "has_prd": plaac.get(header, "0")})
        if len(rows) % 250 == 0:
            print("  %d ..." % len(rows), flush=True)

    out = os.path.join(GYDB, "disorder_all.tsv")
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("wrote %s: %d cores scored; %d skipped as under 20 aa, "
          "%d for non-standard residues" % (out, len(rows), short, nonstandard))


if __name__ == "__main__":
    main()
