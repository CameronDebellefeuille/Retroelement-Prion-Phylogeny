"""Pull the shortlisted Gag candidates out of GyDB and write them up for synthesis.

Nine elements, all PrLD-bearing: six Ty3/Gypsy from insects and three
Retroviridae foamy viruses. Edit CANDIDATES to change the list.

THESE ARE PROTEIN SEQUENCES. GyDB's cores-database is amino acid, and the only
nucleotide data in this repo is the REXdb plant release, which does not contain
any of these elements. So a synthesis order needs one of:

  - codon optimisation from the protein, which the vendor will do (Twist, IDT
    and GenScript all take an amino acid sequence and optimise for a host); or
  - the native nucleotide sequence pulled from GenBank by accession -- GyDB
    lists accessions on its website, they are not in the local files.

Which one matters: codon-optimised DNA gives the same protein but a different
nucleotide sequence, so it is right for expression work and wrong if the
experiment is about the native RNA or the element's own codon usage.

Writes three files to candidates/:

  candidate_gag.faa    full Gag core, one record per element
  candidate_prld.faa   just the PLAAC-called domain
  candidates.tsv       the table, with both sequences and the metadata

Run from the repo root:  python scripts/candidate_selection.py
"""

import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CORES = os.path.join(REPO, "data", "raw", "gydb", "cores-database")
TRAITS = os.path.join(REPO, "data", "processed", "gydb", "traits_272.tsv")
OUT = os.path.join(REPO, "candidates")

CANDIDATES = ["Tom", "297", "17.6", "Yoyo", "HMS-Beagle", "Nomad",
              "EFV", "BFV", "FFV"]


def read_cores(path):
    """GyDB names its cores DOMAIN_element, so GAG_Tom is Tom's Gag."""
    seqs, header = {}, None
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            header = line[1:].strip()
            seqs[header] = []
        elif header and line.strip():
            seqs[header].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    seqs = read_cores(CORES)
    traits = {r["element"]: r for r in
              csv.DictReader(open(TRAITS, encoding="utf-8"), delimiter="\t")}

    missing = [c for c in CANDIDATES if c not in traits]
    if missing:
        sys.exit("not in traits_272.tsv: %s" % ", ".join(missing))

    rows = []
    for name in CANDIDATES:
        t = traits[name]
        gag = seqs["GAG_" + name]
        # PLAAC coordinates are 1-based and inclusive.
        start, end = int(t["prd_start"]), int(t["prd_end"])
        prld = gag[start - 1:end]
        rows.append({
            "element": name,
            "superfamily": t["superfamily"],
            "host": t["host"],
            "gag_len": len(gag),
            "prld_start": start,
            "prld_end": end,
            "prld_len": len(prld),
            "plaac_llr": t["llr"],
            "prld_zone": t["prd_zone"],
            "gag_aa": gag,
            "prld_aa": prld})

    os.makedirs(OUT, exist_ok=True)

    def fasta(path, key, label):
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("; %s -- PROTEIN (amino acid), not DNA. See "
                     "scripts/candidate_selection.py\n" % label)
            for r in rows:
                fh.write(">%s | %s | %s | Gag %d aa | PrLD %d-%d | LLR %s\n%s\n"
                         % (r["element"], r["superfamily"], r["host"] or "host unrecorded",
                            r["gag_len"], r["prld_start"], r["prld_end"],
                            r["plaac_llr"], wrap(r[key])))

    fasta(os.path.join(OUT, "candidate_gag.faa"), "gag_aa", "full Gag cores")
    fasta(os.path.join(OUT, "candidate_prld.faa"), "prld_aa", "PLAAC-called domains only")

    cols = ["element", "superfamily", "host", "gag_len", "prld_start", "prld_end",
            "prld_len", "plaac_llr", "prld_zone", "gag_aa", "prld_aa"]
    with open(os.path.join(OUT, "candidates.tsv"), "w", encoding="utf-8",
              newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print("%d candidates -> %s" % (len(rows), OUT))
    print("%-12s %-13s %-28s %5s %10s %8s" %
          ("element", "superfamily", "host", "Gag", "PrLD", "LLR"))
    for r in rows:
        print("%-12s %-13s %-28s %5d %4d-%-5d %8s"
              % (r["element"], r["superfamily"], r["host"] or "-", r["gag_len"],
                 r["prld_start"], r["prld_end"], r["plaac_llr"]))
    print("\nProtein sequences. For DNA, codon-optimise or pull the native "
          "nucleotide from GenBank -- see the module docstring.")


if __name__ == "__main__":
    main()
