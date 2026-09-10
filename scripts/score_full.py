"""Score the full-length Gag sequences: composition, charge and disorder.

The companion to the PLAAC pass already in data/gag_full/. Same input sequences
-- upstream extension + published core + downstream to the in-frame stop or the
PROT start -- so every measure is on one ruler.

Why not scripts/score.py: that one reads data/gag_*.faa, which is the 75 aa
capped extraction. The cap deletes the elements carrying the signal (capped
PLAAC calls 0 domains in 5,042 copia; uncapped it calls 17), so a disorder or
charge number taken from it is measured on a window rather than on Gag.

Regions are reported separately, using the boundaries gag_full.tsv records:

    [ N-terminus ][ capsid core ][ C-terminal region ]
      gag_upstream  gag_core_len   gag_downstream

The N-terminal columns are the ones the hypothesis is about -- disorder plus
basicity, expected conserved -- so they are kept apart from the whole-sequence
means rather than averaged into them.

Disorder is metapredict's per-residue score; frac_disordered counts residues at
or above 0.5. The GyDB side used the IDR-domain decomposition instead, which
smooths short excursions, so the two frac_disordered columns are close but not
identical and should not be pooled.

metapredict rejects anything outside the 20 standard residues, so sequences
carrying an X get disorder_scored = 0 and empty disorder columns. Composition is
still computed for them, over the non-X residues.

Run under WSL in the copia-prld env, from the repo root:
  python scripts/score_full.py
"""

import csv
import os
from collections import Counter

FULL_DIR = os.path.join("data", "gag_full")
BOUNDS = os.path.join(FULL_DIR, "gag_full.tsv")
PLAAC = os.path.join(FULL_DIR, "gag_full_plaac.tsv")
OUT = os.path.join(FULL_DIR, "traits_full.tsv")

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
BATCH = 500

COLUMNS = ["rexdb_id", "superfamily", "lineage", "species", "is_reference",
           "gag_len", "gag_upstream", "gag_core_len", "gag_downstream",
           "qn", "kr", "de", "net_charge", "gsy",
           "nterm_qn", "nterm_kr", "nterm_de", "nterm_net_charge",
           "disorder_scored", "disorder_mean", "frac_disordered",
           "disorder_nterm", "disorder_core", "disorder_down",
           "has_prd", "llr", "prd_region"]


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
    return {r["rexdb_id"]: r for r in
            csv.DictReader(open(path, encoding="utf-8"), delimiter="\t")}


def composition(seq, prefix=""):
    """Q/N, K/R, D/E, net charge and G/S/Y as fractions of the known residues."""
    known = [c for c in seq if c in STANDARD]
    if not known:
        return {}
    n = len(known)
    c = Counter(known)
    kr = (c["K"] + c["R"]) / n
    de = (c["D"] + c["E"]) / n
    out = {prefix + "qn": (c["Q"] + c["N"]) / n, prefix + "kr": kr,
           prefix + "de": de, prefix + "net_charge": kr - de}
    if not prefix:
        out["gsy"] = (c["G"] + c["S"] + c["Y"]) / n
    return out


def mean(values):
    return sum(values) / len(values) if values else ""


def disorder_rows(seqs):
    """{rexdb_id: per-residue scores} for the sequences metapredict will take."""
    import metapredict as meta

    ids = [i for i in seqs if set(seqs[i]) <= STANDARD]
    print("   disorder: %d of %d scorable (%d carry an X)"
          % (len(ids), len(seqs), len(seqs) - len(ids)), flush=True)
    scores = {}
    for start in range(0, len(ids), BATCH):
        chunk = ids[start:start + BATCH]
        for element, prediction in zip(chunk,
                                       meta.predict_disorder_batch(
                                           [seqs[i] for i in chunk])):
            # batch returns either the score list or (sequence, scores)
            scores[element] = list(
                prediction[1] if isinstance(prediction, (list, tuple))
                and len(prediction) == 2 and not isinstance(prediction[0], float)
                else prediction)
        print("      %d / %d" % (min(start + BATCH, len(ids)), len(ids)),
              flush=True)
    return scores


def main():
    bounds = read_tsv(BOUNDS)
    plaac = read_tsv(PLAAC)
    rows = []

    for short in ("copia", "gypsy"):
        seqs = dict(read_fasta(os.path.join(FULL_DIR,
                                            "gag_full_%s.faa" % short)))
        print("%s: %d sequences" % (short, len(seqs)), flush=True)
        scores = disorder_rows(seqs)

        for element, seq in seqs.items():
            meta_row = bounds[element]
            up = int(meta_row["gag_upstream"] or 0)
            core = int(meta_row["gag_core_len"] or 0)
            hit = plaac.get(element, {})

            row = {c: "" for c in COLUMNS}
            row.update(rexdb_id=element,
                       superfamily=meta_row["superfamily"],
                       lineage=meta_row["lineage"], species=meta_row["species"],
                       is_reference=meta_row["is_reference"],
                       gag_len=len(seq), gag_upstream=up, gag_core_len=core,
                       gag_downstream=meta_row["gag_downstream"],
                       disorder_scored=0,
                       has_prd=hit.get("has_prd", ""), llr=hit.get("llr", ""),
                       prd_region=hit.get("prd_region", ""))
            row.update(composition(seq))
            row.update(composition(seq[:up], "nterm_"))

            d = scores.get(element)
            if d:
                row.update(disorder_scored=1, disorder_mean=mean(d),
                           frac_disordered=sum(1 for x in d if x >= 0.5) / len(d),
                           disorder_nterm=mean(d[:up]),
                           disorder_core=mean(d[up:up + core]),
                           disorder_down=mean(d[up + core:]))
            rows.append(row)

    for row in rows:
        for key, value in row.items():
            if isinstance(value, float):
                row[key] = "%.4f" % value

    with open(OUT, "w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        # reference sequences are named (GEN_Ty1B, YARCTy1-1), the rest are
        # REXdb_ID####; sort the numeric ones numerically and park the names last
        for row in sorted(rows, key=lambda r: (0, int(r["rexdb_id"][9:]), "")
                          if r["rexdb_id"][9:].isdigit()
                          else (1, 0, r["rexdb_id"])):
            writer.writerow(row)
    print("wrote %s (%d rows, %d with disorder)"
          % (OUT, len(rows), sum(1 for r in rows if r["disorder_scored"] == 1)))


if __name__ == "__main__":
    main()
