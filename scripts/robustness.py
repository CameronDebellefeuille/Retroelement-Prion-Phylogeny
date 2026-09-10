"""Rebuild the gypsy RT tree two other ways, to test whether its conclusions are
artefacts of the alignment and trimming settings.

Gypsy is the weaker of the two trees -- 144 columns, 41% of nodes clearing
SH-aLRT 80 / UFBoot 95, and Athila and Ogre recovered at only 46% and 38% purity.
If those numbers move under different settings they were method artefacts; if
they hold, they are properties of gypsy RT.

  baseline   FFT-NS-i + trimAl -automated1   (already built by tree_full.py)
  gappyout   FFT-NS-i + trimAl -gappyout     (trimming sensitivity)
  linsi      MAFFT L-INS-i + -automated1     (alignment sensitivity; L-INS-i is
             what the GyDB tree used and is reachable at 433 tips, though not at
             copia's 688)

Same model set and support settings as tree_full.py, so only the varied stage
differs. Needs the conda environment. Run: python scripts/robustness.py
"""

import os
import subprocess
import time

CLUSTER_DIR = os.path.join("data", "cluster")
OUT_DIR = os.path.join("data", "tree_robust")
REPS = os.path.join(CLUSTER_DIR, "mm_gypsy_rep_seq.fasta")

MODELS = "LG,WAG,JTT,VT,rtREV,HIVb,HIVw,Q.pfam,Dayhoff,Blosum62"

VARIANTS = {
    "gappyout": ("mafft --retree 2 --maxiterate 1000 --thread 8", ["-gappyout"]),
    "linsi": ("mafft --localpair --maxiterate 1000 --thread 8", ["-automated1"]),
}


def run(label, command):
    start = time.time()
    print("   %s ..." % label, flush=True)
    subprocess.run(command, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("   %s done in %.1f min" % (label, (time.time() - start) / 60), flush=True)


def first_seq_len(path):
    length, seen = 0, False
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if seen:
                break
            seen = True
        elif seen:
            length += len(line.strip())
    return length


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for name, (mafft_cmd, trim_flags) in VARIANTS.items():
        prefix = os.path.join(OUT_DIR, name)
        aln, trimmed = prefix + ".aln", prefix + ".trim"
        print("%s:" % name, flush=True)
        run("mafft", ["bash", "-c", "%s '%s' > '%s'" % (mafft_cmd, REPS, aln)])
        run("trimal", ["trimal", "-in", aln, "-out", trimmed] + trim_flags)
        print("   %d columns" % first_seq_len(trimmed), flush=True)
        run("iqtree", ["iqtree", "-s", trimmed, "-m", "MFP", "-mset", MODELS,
                       "-B", "1000", "--alrt", "1000", "-T", "AUTO",
                       "-redo", "--quiet", "-pre", prefix])


if __name__ == "__main__":
    main()
