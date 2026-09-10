"""One RT tree spanning both superfamilies: MAFFT -> occupancy filter -> IQ-TREE.

The separate trees cannot be rooted and cannot be asked whether N-terminal
architecture differs across the copia/gypsy split. This builds the joint tree
that can, on the same 80% cluster representatives -- 688 copia + 433 gypsy.

THE PROBLEM THIS STAGE EXISTS TO SOLVE. REXdb cuts the two superfamilies with
different HMM profiles -- copia RT from Ty1-RT (median 256 aa), gypsy from
Ty3-RT (median 173 aa). Aligned together, the extra copia columns are absent
from every gypsy tip at once, so the missing data is aligned with the deepest
split in the tree and can drive the root by itself.

NEITHER TRIMAL MODE FIXES IT, measured on this alignment (1,121 tips, 466 raw
columns):

    -gappyout      284 columns   gaps: copia 10.5%, gypsy 39.2%
    -automated1    216 columns   gaps: copia 15.9%, gypsy 43.3%

Both leave 80-111 copia-only columns standing, which is the bias itself. A
150+150 probe had suggested trimming removed them cleanly (0.5% / 0.6%); it did
not survive the full set, because both trimAl modes pick their threshold from
the input, so the same rule behaves differently on 300 sequences and on 1,121.
The probe number is not evidence about this alignment and was withdrawn.

THE RULE USED INSTEAD. Keep a column only where at least half of *each*
superfamily has a residue in it. The asymmetry is then removed by construction
rather than by a heuristic that might or might not remove it:

    both-occupancy >= 50%    145 columns   gaps: copia 0.6%, gypsy 0.6%

145 columns is the shared RT core the two profiles genuinely have in common.
The count does not move between a 50% and an 80% threshold -- columns are held
either by nearly all of both superfamilies or by one alone -- so the rule is not
balanced on a knife edge. For scale, the GyDB tree recovered four superfamilies
as contiguous blocks on 131 columns.

What is lost is resolution, not correctness: copia alone supports 185 columns.
Read clade membership from this tree, and the deep backbone only where support
allows.

Rooting stays an open decision. This tree makes rooting *possible*; it does not
assert one. REXdb's own outgroup lineages are present as tips in both sets.

Long running. Needs the conda environment. Run from the repo root:
  python scripts/tree_joint.py
"""

import csv
import os
import subprocess
import time

CLUSTER = os.path.join("data", "cluster")
TREE_DIR = os.path.join("data", "tree_joint")
PREFIX = os.path.join(TREE_DIR, "joint")
ELEMENTS = os.path.join("data", "elements.tsv")

MIN_OCCUPANCY = 0.5     # of each superfamily, per column
MODELS = "LG,WAG,JTT,VT,rtREV,HIVb,HIVw,Q.pfam,Dayhoff,Blosum62"
UFBOOT = 1000
ALRT = 1000


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


def run(label, command):
    start = time.time()
    print("   %s ..." % label, flush=True)
    subprocess.run(command, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("   %s done in %.1f min" % (label, (time.time() - start) / 60), flush=True)


def occupancy_filter(aln, out, superfamily):
    """Keep columns where each superfamily is at least MIN_OCCUPANCY occupied."""
    records = list(read_fasta(aln))
    width = len(records[0][1])
    groups = {}
    for element, seq in records:
        groups.setdefault(superfamily[element], []).append(seq)

    keep = []
    for i in range(width):
        if all(sum(1 for s in seqs if s[i] != "-") / len(seqs) >= MIN_OCCUPANCY
               for seqs in groups.values()):
            keep.append(i)

    with open(out, "w", encoding="utf-8", newline="\n") as handle:
        for element, seq in records:
            handle.write(">%s\n%s\n" % (element, "".join(seq[i] for i in keep)))

    for name, seqs in sorted(groups.items()):
        gaps = sum(sum(1 for i in keep if s[i] == "-") for s in seqs)
        print("   %-10s %4d tips, %.1f%% gaps in the kept columns"
              % (name, len(seqs), 100 * gaps / (len(seqs) * len(keep))), flush=True)
    return width, len(keep)


def main():
    os.makedirs(TREE_DIR, exist_ok=True)
    superfamily = {r["rexdb_id"]: r["superfamily"] for r in
                   csv.DictReader(open(ELEMENTS, encoding="utf-8"),
                                  delimiter="\t")}

    combined = PREFIX + ".faa"
    counts, seen = {}, set()
    with open(combined, "w", encoding="utf-8", newline="\n") as handle:
        for short in ("copia", "gypsy"):
            reps = os.path.join(CLUSTER, "mm_%s_rep_seq.fasta" % short)
            counts[short] = 0
            for element, seq in read_fasta(reps):
                # REXdb IDs are unique across the whole release, so the two sets
                # can share one file without prefixing -- which keeps the tip
                # labels joinable to elements.tsv and traits_full.tsv unchanged
                if element in seen:
                    raise SystemExit("duplicate tip id across superfamilies: %s"
                                     % element)
                seen.add(element)
                handle.write(">%s\n%s\n" % (element, seq))
                counts[short] += 1
    print("%d tips: copia %d, gypsy %d"
          % (sum(counts.values()), counts["copia"], counts["gypsy"]), flush=True)

    aln, trimmed = PREFIX + ".aln", PREFIX + ".trim"
    run("mafft", ["bash", "-c",
                  "mafft --retree 2 --maxiterate 1000 --thread 8 '%s' > '%s'"
                  % (combined, aln)])

    print("   occupancy filter ...", flush=True)
    before, after = occupancy_filter(aln, trimmed, superfamily)
    print("   alignment %d -> %d columns" % (before, after), flush=True)

    run("iqtree", ["iqtree", "-s", trimmed, "-m", "MFP", "-mset", MODELS,
                   "-B", str(UFBOOT), "--alrt", str(ALRT),
                   "-T", "AUTO", "-redo", "--quiet", "-pre", PREFIX])
    print("tree at %s.treefile" % PREFIX, flush=True)


if __name__ == "__main__":
    main()
