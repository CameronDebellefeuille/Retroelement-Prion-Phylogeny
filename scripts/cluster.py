"""Cluster REXdb RT at 80% identity to pick tree tips, the way the source paper did.

Neumann et al. 2019 (Mobile DNA 10:1, doi 10.1186/s13100-018-0144-1) reduced
redundancy by clustering sequences that "shared at least 80% identity over at
least 90% of their length", then took one representative per cluster. That is
reproduced here with MMseqs2 rather than the paper's CD-HIT -- MMseqs2 is already
in environment.yml and CD-HIT is not. The threshold pair is the paper's.

Copia and gypsy are clustered separately and get separate trees. Gypsy RT is cut
by a different HMM profile (median 173 aa vs 256), so a joint alignment puts an
~83-column gap block in every gypsy tip; see F-2 and the note in tree.py.

Filters, all fixed before running per L-8:

  rt_ambiguous == 0   X in RT corrupts both the clustering identity and the
                      alignment columns downstream. Drops ~20% copia, ~16% gypsy.
  has_gag == 1        every tip must be able to carry the PrLD trait, which is
  gag_ambiguous == 0  the GyDB design -- its 272-element set required a clean
                      GAG and a clean RT for the same element.

Also joins the full REXdb lineage path, which elements.tsv does not carry:
that file stops at column 5 of the classification, so all 8,749 gypsy elements
collapse to "chromovirus"/"non-chromovirus" and the real lineage names (Retand,
Tekay, Athila, CRM, Reina, Ogre, ...) are lost. They live in columns 6-8.

Needs the conda environment. Run: python scripts/cluster.py
"""

import csv
import os
import subprocess
from collections import defaultdict

OUT = "data"
CLUSTER_DIR = os.path.join(OUT, "cluster")
CLASSIFICATION = os.path.join(OUT, "REXdb", "Viridiplantae_v4.0.classification")

MIN_SEQ_ID = 0.8      # paper: "at least 80% identity"
MIN_COVERAGE = 0.9    # paper: "over at least 90% of their length"
COV_MODE = 0          # bidirectional, i.e. 90% of both sequences
THREADS = 8

SUPERFAMILIES = {"Ty1/copia": "copia", "Ty3/gypsy": "gypsy"}


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


def read_lineages():
    """rexdb_id -> (full path, leaf lineage) for LTR elements."""
    lineages = {}
    for row in open(CLASSIFICATION, encoding="utf-8"):
        fields = row.rstrip("\n").split("\t")
        if len(fields) < 5 or fields[2] != "LTR":
            continue
        path = "/".join(fields[3:])
        lineages[fields[0]] = (path, fields[-1])
    return lineages


def keep(row):
    return (row["rt_ambiguous"] == "0"
            and row["has_gag"] == "1"
            and row["gag_ambiguous"] == "0")


def cluster(short, source, ids):
    """Write the filtered set, cluster it, return representative -> members."""
    filtered = os.path.join(CLUSTER_DIR, "rt_%s_filtered.faa" % short)
    with open(filtered, "w", encoding="utf-8", newline="\n") as handle:
        for element, seq in read_fasta(source):
            if element in ids:
                # B/J break trimAl's similarity matrix; X is the standard unknown
                handle.write(">%s\n%s\n"
                             % (element, seq.replace("B", "X").replace("J", "X")))

    prefix = os.path.join(CLUSTER_DIR, "mm_%s" % short)
    tmp = os.path.join(CLUSTER_DIR, "tmp_%s" % short)
    print("   mmseqs %s ..." % short, flush=True)
    subprocess.run(["mmseqs", "easy-cluster", filtered, prefix, tmp,
                    "--min-seq-id", str(MIN_SEQ_ID),
                    "-c", str(MIN_COVERAGE),
                    "--cov-mode", str(COV_MODE),
                    "--threads", str(THREADS), "-v", "0"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    members = defaultdict(list)
    with open(prefix + "_cluster.tsv", encoding="utf-8") as handle:
        for line in handle:
            rep, member = line.rstrip("\n").split("\t")
            members[rep].append(member)
    return members


def main():
    os.makedirs(CLUSTER_DIR, exist_ok=True)
    elements = {r["rexdb_id"]: r for r in
                csv.DictReader(open(os.path.join(OUT, "elements.tsv"),
                                    encoding="utf-8"), delimiter="\t")}
    lineages = read_lineages()

    summary = []
    for superfamily, short in SUPERFAMILIES.items():
        pool = [r for r in elements.values() if r["superfamily"] == superfamily]
        usable = {r["rexdb_id"] for r in pool if keep(r)}
        print("%s: %d elements, %d pass filters"
              % (superfamily, len(pool), len(usable)))

        members = cluster(short, os.path.join(OUT, "rt_%s.faa" % short), usable)
        print("   %d clusters at %d%% id / %d%% coverage"
              % (len(members), MIN_SEQ_ID * 100, MIN_COVERAGE * 100))

        tips = os.path.join(CLUSTER_DIR, "tips_%s.tsv" % short)
        with open(tips, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["rexdb_id", "superfamily", "lineage_elements",
                             "lineage_leaf", "lineage_path", "species", "taxid",
                             "rt_len", "gag_len", "gag_upstream", "cluster_size",
                             "members"])
            for rep in sorted(members, key=lambda r: -len(members[r])):
                row = elements[rep]
                path, leaf = lineages.get(rep, ("NA", "NA"))
                writer.writerow([rep, superfamily, row["lineage"], leaf, path,
                                 row["species"], row["taxid"], row["rt_len"],
                                 row["gag_len"], row["gag_upstream"],
                                 len(members[rep]), ",".join(sorted(members[rep]))])

        # count tips per true lineage, which is the resolution elements.tsv lacks
        per_leaf = defaultdict(int)
        for rep in members:
            per_leaf[lineages.get(rep, ("NA", "NA"))[1]] += 1
        print("   lineages: %s" % ", ".join(
            "%s %d" % kv for kv in sorted(per_leaf.items(), key=lambda k: -k[1])))
        summary.append((superfamily, len(pool), len(usable), len(members)))

    print()
    for superfamily, total, usable, clusters in summary:
        print("%-12s %5d elements -> %5d after filters -> %4d tips"
              % (superfamily, total, usable, clusters))


if __name__ == "__main__":
    main()
