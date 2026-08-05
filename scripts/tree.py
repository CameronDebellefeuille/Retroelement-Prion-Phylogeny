"""Build the copia RT phylogeny: sample -> MAFFT -> trimAl -> IQ-TREE.

Sampling (L-2) keeps every element with a long N-terminal extension, then fills
out the backbone evenly across lineages. The long tail is where the trait signal
lives -- 174 copia elements carry >=180 aa upstream and they hold nearly all the
prion-like domains -- so balanced sampling alone would drop the signal.

Copia only. Gypsy RT is cut by a different HMM profile (median 173 aa vs 256),
so a joint alignment would put an ~83-column gap block in every gypsy tip; see
F-2. Rooting is deferred (L-3).

Fast settings: LG+G4, SH-aLRT support, no ModelFinder, no bootstraps.

Needs the conda environment. Run: python scripts/tree.py
"""

import csv
import os
import random
import subprocess
from collections import defaultdict

OUT = "data"
TREE_DIR = os.path.join(OUT, "tree")
SUPERFAMILY = "Ty1/copia"

KEEP_ALL_UPSTREAM = 180   # every element at least this long is retained
TARGET_TIPS = 700
SEED = 1


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


def choose_tips(elements):
    """Every long-N-terminus element, then an even fill across lineages."""
    usable = [r for r in elements.values()
              if r["superfamily"] == SUPERFAMILY and r["in_tree"] == "1"]
    long_tail = {r["rexdb_id"] for r in usable
                 if r["gag_upstream"] and int(r["gag_upstream"]) >= KEEP_ALL_UPSTREAM}

    by_lineage = defaultdict(list)
    for r in usable:
        if r["rexdb_id"] not in long_tail:
            by_lineage[r["lineage"]].append(r["rexdb_id"])

    random.seed(SEED)
    for ids in by_lineage.values():
        ids.sort()
        random.shuffle(ids)

    # round-robin across lineages so small ones are not swamped
    chosen = set(long_tail)
    order = sorted(by_lineage, key=lambda k: len(by_lineage[k]))
    i = 0
    while len(chosen) < TARGET_TIPS and any(by_lineage.values()):
        pool = by_lineage[order[i % len(order)]]
        if pool:
            chosen.add(pool.pop())
        i += 1
    return chosen, long_tail


def run(label, command):
    print("   %s ..." % label, flush=True)
    subprocess.run(command, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main():
    os.makedirs(TREE_DIR, exist_ok=True)
    elements = {r["rexdb_id"]: r for r in
                csv.DictReader(open(os.path.join(OUT, "elements.tsv"),
                                    encoding="utf-8"), delimiter="\t")}
    chosen, long_tail = choose_tips(elements)
    print("sampled %d tips (%d with >=%d aa upstream, kept unconditionally)"
          % (len(chosen), len(long_tail), KEEP_ALL_UPSTREAM))

    lineages = defaultdict(int)
    for i in chosen:
        lineages[elements[i]["lineage"]] += 1
    print("   lineages: %s" % ", ".join("%s %d" % kv for kv in
                                        sorted(lineages.items(), key=lambda k: -k[1])))

    sampled = os.path.join(TREE_DIR, "rt_sampled.faa")
    with open(sampled, "w", encoding="utf-8", newline="\n") as handle:
        for element, seq in read_fasta(os.path.join(OUT, "rt_copia.faa")):
            if element in chosen:
                # B/J break trimAl's similarity matrix; X is the standard unknown
                handle.write(">%s\n%s\n" % (element, seq.replace("B", "X").replace("J", "X")))

    aln = os.path.join(TREE_DIR, "rt.aln")
    trimmed = os.path.join(TREE_DIR, "rt.trim")
    run("mafft", ["bash", "-c", "mafft --auto --thread 4 %s > %s" % (sampled, aln)])
    run("trimal", ["trimal", "-in", aln, "-out", trimmed, "-gappyout"])
    run("iqtree", ["iqtree", "-s", trimmed, "-m", "LG+G4", "-T", "AUTO",
                   "--alrt", "1000", "-redo", "--quiet",
                   "-pre", os.path.join(TREE_DIR, "rt")])

    columns = 0
    for header, seq in read_fasta(trimmed):
        columns = len(seq)
        break
    print("alignment %d columns; tree at %s.treefile" % (columns, os.path.join(TREE_DIR, "rt")))


if __name__ == "__main__":
    main()
