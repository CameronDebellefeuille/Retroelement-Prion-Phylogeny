"""Build the full REXdb RT phylogenies from the 80% clusters: MAFFT -> trimAl -> IQ-TREE.

One tree per superfamily, on the cluster representatives from cluster.py. Copia
and gypsy are never aligned together: gypsy RT is cut by a different HMM profile
(median 173 aa vs 256), which would put an ~83-column gap block in every gypsy
tip. See F-2 and tree.py.

Settings follow the GyDB tree so the two are comparable, with one deliberate
deviation:

  align    MAFFT FFT-NS-i (--retree 2 --maxiterate 1000). The GyDB tree used
           L-INS-i, whose practical ceiling is a few hundred sequences; at 688
           copia tips in 7 GB that is not reachable. Deviation is scale-forced.
  trim     trimAl -automated1, as GyDB.
  model    IQ-TREE ModelFinder over the same curated candidate set as GyDB. The
           full sweep makes no progress in reasonable time; rtREV is included as
           the domain-appropriate candidate estimated from retroviral RT.
  support  1000 UFBoot + 1000 SH-aLRT, as GyDB.

Rooting is deferred (L-3). REXdb carries its own outgroup lineages
(Ty1-outgroup, chromo-outgroup, non-chromo-outgroup) which are present as tips.

Long running -- IQ-TREE on 688 tips is the bottleneck. Needs the conda
environment. Run: python scripts/tree_full.py [copia|gypsy]
"""

import os
import subprocess
import sys
import time

OUT = "data"
CLUSTER_DIR = os.path.join(OUT, "cluster")
TREE_DIR = os.path.join(OUT, "tree_full")

MODELS = "LG,WAG,JTT,VT,rtREV,HIVb,HIVw,Q.pfam,Dayhoff,Blosum62"
UFBOOT = 1000
ALRT = 1000
THREADS = "AUTO"


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


def build(short):
    reps = os.path.join(CLUSTER_DIR, "mm_%s_rep_seq.fasta" % short)
    tips = sum(1 for _ in read_fasta(reps))
    prefix = os.path.join(TREE_DIR, short)
    aln = prefix + ".aln"
    trimmed = prefix + ".trim"

    print("%s: %d tips" % (short, tips), flush=True)
    run("mafft", ["bash", "-c",
                  "mafft --retree 2 --maxiterate 1000 --thread 8 '%s' > '%s'"
                  % (reps, aln)])

    before = len(next(read_fasta(aln))[1])
    run("trimal", ["trimal", "-in", aln, "-out", trimmed, "-automated1"])
    after = len(next(read_fasta(trimmed))[1])
    print("   alignment %d -> %d columns" % (before, after), flush=True)

    run("iqtree", ["iqtree", "-s", trimmed, "-m", "MFP", "-mset", MODELS,
                   "-B", str(UFBOOT), "--alrt", str(ALRT),
                   "-T", THREADS, "-redo", "--quiet", "-pre", prefix])
    print("   tree at %s.treefile" % prefix, flush=True)


def main():
    os.makedirs(TREE_DIR, exist_ok=True)
    targets = sys.argv[1:] or ["gypsy", "copia"]   # smaller first, so failures surface early
    for short in targets:
        build(short)


if __name__ == "__main__":
    main()
