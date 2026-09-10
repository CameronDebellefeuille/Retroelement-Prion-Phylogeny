"""Identify the nucleocapsid zinc knuckle in the full-length Gag sequences, by HMM.

Replaces the regex `C-X2-C-X4-H-X4-C` used in earlier probes. That pattern is
essentially the classic PROSITE-style consensus, and on this data it calls a
knuckle in 70% of copia Gags where the Pfam models call one in 21% -- 2,520
regex-only calls against 6 Pfam-only. A 14-residue pattern with four fixed
positions hits by chance in a 300 aa protein; a profile does not. Everything
CCHC-based should use this table, not the regex.

Sensitivity was the worry, so the whole Pfam zinc-knuckle clan is used rather
than PF00098 alone, and each model is run twice:

    strict      --cut_ga, Pfam's own curated gathering threshold
    sensitive   -E 1e-3, to catch knuckles too divergent for the threshold

Both counts are written per element, so any analysis can state which definition
it used and check that its result does not depend on the choice.

Not done here, and worth knowing: a full InterProScan would add PROSITE, SMART,
CDD and PRINTS evidence for the same domain. The web service cannot take 13,572
sequences, and a local install is a large download. The clan sweep plus the
permissive threshold is the cheap approximation to it; if a knuckle call ever
becomes load-bearing on its own, run InterProScan on that subset.

Needs the conda environment (hmmer). Run from the repo root:

  python scripts/nc_scan.py                          # the REXdb Gag sets
  python scripts/nc_scan.py <fasta> <out.tsv>        # any other protein set

The second form is how the GyDB cores are scanned, so both databases get their
knuckle calls from the same models at the same thresholds.
"""

import gzip
import os
import subprocess
import sys
import urllib.request
from collections import defaultdict

FULL = os.path.join("data", "gag_full")
HMM_DIR = os.path.join(FULL, "pfam")
PANEL = os.path.join(HMM_DIR, "nc_panel.hmm")
OUT = os.path.join(FULL, "nc_hits.tsv")

# Pfam clan CL0511, restricted to the zinc-knuckle families. The clan also holds
# eIF3G, CLIP1, AIR2 and Lin-28A knuckles, which are host proteins; they are left
# out so a hit here means a retroelement-style knuckle.
MODELS = {
    "PF00098": "zf-CCHC",
    "PF13696": "zf-CCHC_2",
    "PF13917": "zf-CCHC_3",
    "PF14392": "zf-CCHC_4",
    "PF14787": "GAG_viral_zf",
    "PF15288": "zf-CCHC_6",
}

SENSITIVE_E = "1e-3"

# PROSITE PS50158 run locally with pftools. Validated against InterProScan on the
# 336 GyDB cores: Pfam(E<1e-3) alone recovers 87% of InterPro's knuckle calls,
# PROSITE alone 76%, and the union of the two 92% at 99% specificity. That union
# is the closest thing to InterPro that can be run over 13,572 sequences, and is
# what `nc_combined` records.
PROSITE_PRF = os.path.expanduser("~/prosite/ps50158.prf")


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


def build_panel():
    os.makedirs(HMM_DIR, exist_ok=True)
    if os.path.exists(PANEL):
        return
    chunks = []
    for accession, name in MODELS.items():
        path = os.path.join(HMM_DIR, accession + ".hmm")
        if not os.path.exists(path):
            url = ("https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/%s"
                   "?annotation=hmm" % accession)
            with urllib.request.urlopen(url, timeout=60) as response:
                text = gzip.decompress(response.read()).decode()
            open(path, "w", encoding="utf-8", newline="\n").write(text)
            print("   fetched %s %s" % (accession, name), flush=True)
        chunks.append(open(path, encoding="utf-8").read())
    open(PANEL, "w", encoding="utf-8", newline="\n").write("".join(chunks))
    subprocess.run(["hmmpress", "-f", PANEL], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def scan(fasta, tag, threshold):
    out = os.path.join(HMM_DIR, "%s_%s.domtbl" % (tag, threshold))
    flag = ["--cut_ga"] if threshold == "strict" else ["-E", SENSITIVE_E]
    subprocess.run(["hmmsearch"] + flag + ["--domtblout", out, "--cpu", "4",
                                           PANEL, fasta],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    hits = defaultdict(list)
    for line in open(out, encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.split()
        # target, model accession, envelope start/end on the target, i-Evalue
        hits[f[0]].append((f[4].split(".")[0], int(f[19]), int(f[20]), float(f[12])))
    return hits


def prosite(fasta):
    """Sequences with a PS50158 profile match, or None if pftools is absent."""
    if not os.path.exists(PROSITE_PRF):
        return None
    try:
        out = subprocess.run(["pfsearch", "-f", PROSITE_PRF, fasta],
                             capture_output=True, text=True, check=True).stdout
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    hits = set()
    for line in out.splitlines():
        f = line.split()
        if len(f) >= 7 and f[2] == "pos.":
            hits.add(f[6])
    return hits


def main():
    print("building the knuckle panel ...", flush=True)
    build_panel()

    if len(sys.argv) == 3:
        sets, out_path = [(sys.argv[1], os.path.basename(sys.argv[1]))], sys.argv[2]
    else:
        sets = [(os.path.join(FULL, "gag_full_%s.faa" % s), s)
                for s in ("copia", "gypsy")]
        out_path = OUT

    rows = {}
    for fasta, short in sets:
        ids = [element for element, _ in read_fasta(fasta)]
        print("scanning %s (%d sequences) ..." % (short, len(ids)), flush=True)
        strict = scan(fasta, short, "strict")
        loose = scan(fasta, short, "sensitive")
        pro = prosite(fasta)
        if pro is None:
            print("   PROSITE skipped (pftools or the profile is not installed)",
                  flush=True)
        for element in ids:
            s, l = strict.get(element, []), loose.get(element, [])
            best = min(s or l, key=lambda h: h[3]) if (s or l) else None
            rows[element] = {
                "rexdb_id": element,
                "nc_strict": int(bool(s)),
                "nc_sensitive": int(bool(l)),
                "n_strict": len(s),
                "n_sensitive": len(l),
                "best_model": best[0] if best else "",
                "best_start": best[1] if best else "",
                "best_end": best[2] if best else "",
                "best_evalue": "%.2g" % best[3] if best else "",
                "nc_prosite": "" if pro is None else int(element in pro),
                "nc_combined": "" if pro is None
                               else int(bool(l) or element in pro),
            }
        print("   strict %d, sensitive %d, PROSITE %s, combined %s -- of %d"
              % (sum(1 for e in ids if rows[e]["nc_strict"]),
                 sum(1 for e in ids if rows[e]["nc_sensitive"]),
                 "-" if pro is None else sum(1 for e in ids if rows[e]["nc_prosite"]),
                 "-" if pro is None else sum(1 for e in ids if rows[e]["nc_combined"]),
                 len(ids)), flush=True)

    columns = ["rexdb_id", "nc_strict", "nc_sensitive", "nc_prosite",
               "nc_combined", "n_strict", "n_sensitive",
               "best_model", "best_start", "best_end", "best_evalue"]
    with open(out_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(columns) + "\n")
        for element in sorted(rows, key=lambda e: (len(e), e)):
            handle.write("\t".join(str(rows[element][c]) for c in columns) + "\n")
    print("\nwrote %s (%d rows)" % (out_path, len(rows)))


if __name__ == "__main__":
    main()
