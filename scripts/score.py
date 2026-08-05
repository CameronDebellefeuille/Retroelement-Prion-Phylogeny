"""Score every Gag sequence for the trait measures, into one table.

Composition and charge are computed here. PLAAC is called as a jar (alpha 1.0 is
primary, 0 and 0.5 recorded for sensitivity -- L-4). Disorder comes from
metapredict (L-6), which rejects X, so sequences carrying an unresolved residue
get an empty disorder column rather than a substituted one.

gag_upstream marks where the capsid core begins, so disorder is reported for the
N-terminal region and the core separately.

Needs the conda environment (java, metapredict). Run: python scripts/score.py
"""

import csv
import os
import subprocess
from collections import Counter

OUT = "data"
PLAAC_JAR = os.path.join("plaac", "plaac.jar")
CORE_LENGTH = 60          # L-5
ALPHAS = ["1.0", "0.5", "0"]  # first is primary

COLUMNS = ["rexdb_id", "superfamily", "lineage", "species", "is_reference",
           "gag_len", "gag_upstream", "qn", "kr", "de", "net_charge", "gsy",
           "llr", "nllr", "prd_score", "prd_start", "prd_end", "has_prd",
           "llr_a05", "llr_a0",
           "disorder_scored", "disorder_mean", "disorder_nterm",
           "disorder_core", "disorder_fraction"]


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


def composition(seq):
    n = len(seq)
    c = Counter(seq)
    kr = (c["K"] + c["R"]) / n
    de = (c["D"] + c["E"]) / n
    return {"qn": (c["Q"] + c["N"]) / n, "kr": kr, "de": de,
            "net_charge": kr - de, "gsy": (c["G"] + c["S"] + c["Y"]) / n}


def run_plaac(fasta, alpha):
    """{rexdb_id: {llr, nllr, prd_score, prd_start, prd_end, has_prd}}"""
    result = subprocess.run(
        ["java", "-Xmx4g", "-jar", PLAAC_JAR, "-i", fasta,
         "-a", alpha, "-c", str(CORE_LENGTH)],
        capture_output=True, text=True, check=True)
    lines = [l for l in result.stdout.splitlines() if not l.startswith("#")]
    header = lines[0].split("\t")
    at = {c: header.index(c) for c in
          ("SEQid", "LLR", "NLLR", "PRDscore", "PRDstart", "PRDend", "COREscore")}
    out = {}
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) < len(header):
            continue
        called = f[at["COREscore"]] not in ("NaN", "")
        out[f[at["SEQid"]]] = {
            "llr": f[at["LLR"]], "nllr": f[at["NLLR"]],
            "prd_score": f[at["PRDscore"]] if called else "",
            "prd_start": f[at["PRDstart"]] if called else "",
            "prd_end": f[at["PRDend"]] if called else "",
            "has_prd": int(called)}
    return out


def run_metapredict(seqs, upstream):
    """{rexdb_id: {...}} for sequences without X. metapredict rejects X."""
    import metapredict as meta
    ids = [i for i in seqs if "X" not in seqs[i]]
    scores = meta.predict_disorder_batch([seqs[i] for i in ids])
    out = {}
    for element, prediction in zip(ids, scores):
        # batch returns either the score list or (sequence, scores)
        d = list(prediction[1] if isinstance(prediction, (list, tuple))
                 and len(prediction) == 2 and not isinstance(prediction[0], float)
                 else prediction)
        up = upstream.get(element, 0)
        out[element] = {
            "disorder_scored": 1,
            "disorder_mean": sum(d) / len(d),
            "disorder_fraction": sum(1 for x in d if x >= 0.5) / len(d),
            "disorder_nterm": sum(d[:up]) / up if 0 < up < len(d) else "",
            "disorder_core": sum(d[up:]) / (len(d) - up) if 0 < up < len(d) else ""}
    return out


def main():
    elements = {r["rexdb_id"]: r for r in
                csv.DictReader(open(os.path.join(OUT, "elements.tsv"),
                                    encoding="utf-8"), delimiter="\t")}
    rows = []
    for short in ("copia", "gypsy"):
        fasta = os.path.join(OUT, "gag_%s.faa" % short)
        seqs = dict(read_fasta(fasta))
        upstream = {i: int(elements[i]["gag_upstream"] or 0) for i in seqs}
        print("%s: %d sequences" % (short, len(seqs)))

        plaac = {a: run_plaac(fasta, a) for a in ALPHAS}
        print("   plaac done")
        disorder = run_metapredict(seqs, upstream)
        print("   disorder done (%d scored, %d skipped for X)"
              % (len(disorder), len(seqs) - len(disorder)))

        for element, seq in seqs.items():
            meta_row = elements[element]
            row = {c: "" for c in COLUMNS}
            row.update(rexdb_id=element, superfamily=meta_row["superfamily"],
                       lineage=meta_row["lineage"], species=meta_row["species"],
                       is_reference=meta_row["is_reference"],
                       gag_len=len(seq), gag_upstream=upstream[element],
                       disorder_scored=0)
            row.update(composition(seq))
            row.update(plaac[ALPHAS[0]].get(element, {}))
            row["llr_a05"] = plaac["0.5"].get(element, {}).get("llr", "")
            row["llr_a0"] = plaac["0"].get(element, {}).get("llr", "")
            row.update(disorder.get(element, {}))
            rows.append(row)

    path = os.path.join(OUT, "traits.tsv")
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in sorted(rows, key=lambda r: int(r["rexdb_id"].replace("REXdb_ID", ""))):
            writer.writerow(row)
    print("wrote %s (%d rows)" % (path, len(rows)))


if __name__ == "__main__":
    main()
