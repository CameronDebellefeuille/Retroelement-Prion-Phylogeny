"""metapredict disorder for the 272 Gag cores, to map beside the PLAAC calls.

Uses metapredict's own IDR decomposition (predict_disorder_domains) rather than
a hand-picked score cutoff -- the tool's thresholds are calibrated, an invented
one would be another undecided parameter, and L-6 only settled the tool.

Writes disorder_272.tsv: mean score, fraction of the core in an IDR, longest IDR,
and the mean score over the first 40 residues, which is the N-terminal end the
project's hypothesis is actually about.

Run under WSL in the copia-prld env:
  python scripts/disorder_272.py
"""

import csv
import os

import metapredict as meta

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
GYDB = os.path.join(REPO, "data", "processed", "gydb")
RAW  = os.path.join(REPO, "data", "raw", "gydb")


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


def main():
    rows = []
    for name, seq in read_fasta(os.path.join(GYDB, "gag_272.faa")):
        result = meta.predict_disorder_domains(seq)
        scores = result.disorder
        idrs = result.disordered_domain_boundaries
        in_idr = sum(end - start for start, end in idrs)
        longest = max((end - start for start, end in idrs), default=0)
        head = scores[:40]
        rows.append({
            "tip": name,
            "gag_len": len(seq),
            "mean_disorder": round(sum(scores) / len(scores), 4),
            "frac_disordered": round(in_idr / len(seq), 4),
            "longest_idr": longest,
            "n_idr": len(idrs),
            "nterm40_disorder": round(sum(head) / len(head), 4)})
        if len(rows) % 50 == 0:
            print("  %d ..." % len(rows), flush=True)

    out = os.path.join(GYDB, "disorder_272.tsv")
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    frac = sorted(r["frac_disordered"] for r in rows)
    print("wrote %s for %d sequences" % (out, len(rows)))
    print("fraction disordered: median %.2f, IQR %.2f-%.2f, min %.2f, max %.2f"
          % (frac[len(frac) // 2], frac[len(frac) // 4], frac[3 * len(frac) // 4],
             frac[0], frac[-1]))
    print("entirely ordered (no IDR at all): %d"
          % sum(1 for r in rows if r["n_idr"] == 0))


if __name__ == "__main__":
    main()
