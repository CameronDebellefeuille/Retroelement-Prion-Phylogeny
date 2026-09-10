"""Assemble the 272-element GyDB set: RT for the tree, Gag traits for the tips.

The set is every element in the cores database carrying both a GAG and an RT
core with no X in either. Same design principle as the pipeline -- RT is
aligned to build the phylogeny, Gag is scored to make the trait, and the two are
joined by element name. The tree is never built from Gag.

PLAAC values come from the whole-database alpha=0.5 run (all_cores_a05.tsv), not
from a 272-sequence rerun. At alpha < 1 the background is partly the input's own
composition, so scoring the subset alone would give different numbers from every
figure reported so far; scoring once over all 2,636 cores keeps one ruler.

Writes rt_272.faa (tree input), gag_272.faa, and traits_272.tsv (tip data).
Run from the repo root: python scripts/prep_272.py
"""

import collections
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
GYDB = os.path.join(REPO, "data", "processed", "gydb")
RAW  = os.path.join(REPO, "data", "raw", "gydb")
CORES = os.path.join(RAW, "cores-database")
# Cameron's curated classification sheet. It lived in ~/Downloads until it was
# moved into the repo -- it is the only source of superfamily, host and clade on
# the tree, there is no other copy, and it is not a GyDB download. Treat it as
# raw data, not as a derived file.
SUPERFAMILY_CSV = os.path.join(RAW, "gydb_superfamily_host.csv")
TWO = {"PR_ULP1", "PR_OTU", "ORF1_Nterdomain"}   # two-token abbreviations
SHORT_RT = 150   # below this an RT core is mostly gap in the alignment


def read_cores(path):
    seqs, header = {}, None
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            header = line[1:].strip()
            seqs[header] = []
        elif header and line.strip():
            seqs[header].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def split_header(header):
    """GAG_Ty1B -> (GAG, Ty1B); PR_ULP1_HmGINA1 -> (PR_ULP1, HmGINA1)."""
    parts = header.split("_")
    two = "_".join(parts[:2])
    if two in TWO:
        return two, "_".join(parts[2:])
    return parts[0], "_".join(parts[1:])


def tip_label(name):
    """Newick-safe: CopiaSL_monotypic|Chr08_2s54 -> CopiaSL_monotypic-Chr08_2s54."""
    return re.sub(r"[^A-Za-z0-9._-]", "-", name)


def group_of(name):
    if "monotypic" in name:
        return "monotypic"
    if "CopiaSL_" in name or "GypsySL_" in name:
        return "superlineage"
    return "named"


def read_superfamilies(path):
    """Cameron's GyDB classification sheet: element name -> superfamily, host, clade.

    169 of its 437 rows have no superfamily -- the GIN/GINGER/BIN/CIN entries,
    which are not LTR elements -- and are skipped. FFV is filed as "Retrovirus"
    rather than "Retroviridae"; it is a spumaretrovirus, so the two are merged.
    """
    table = {}
    for row in csv.DictReader(open(path, encoding="utf-8-sig")):
        family = row["Superfamily/Family"].strip()
        if not family:
            continue
        if family == "Retrovirus":
            family = "Retroviridae"
        table[row["Element Name"].strip()] = {
            "superfamily": family,
            "host": row["Host Species"].strip(),
            "clade": row["Clade"].strip()}
    return table


def classify(name, table, lowered):
    """Match a tip to the sheet, then fall back to what the name itself says."""
    if name in table:
        return table[name]
    if name.lower() in lowered:                       # GALV vs GaLV
        return table[lowered[name.lower()]]
    stem = re.sub(r"[\d.-]+$", "", name)              # Calypso5-1 vs Calypso
    if stem and stem in table:
        return table[stem]
    if name.startswith("CopiaSL"):                    # the name is the assignment
        return {"superfamily": "Ty1/Copia", "host": "", "clade": ""}
    if name.startswith("GypsySL"):
        return {"superfamily": "Ty3/Gypsy", "host": "", "clade": ""}
    return {"superfamily": "unclassified", "host": "", "clade": ""}


CCHC = re.compile(r"C.{2}C.{4}H.{4}C")   # retroviral zinc knuckle, F-2's pattern


def prd_zone(mid_pct):
    """Which third of the protein the called domain sits in."""
    if mid_pct is None:
        return "no PrLD"
    if mid_pct < 100 / 3:
        return "N-terminal"
    if mid_pct > 200 / 3:
        return "C-terminal"
    return "middle"


def main():
    seqs = read_cores(CORES)
    elements = collections.defaultdict(dict)
    for header in seqs:
        protein, name = split_header(header)
        elements[name][protein] = header

    keep = sorted(name for name, domains in elements.items()
                  if "GAG" in domains and "RT" in domains
                  and "X" not in seqs[domains["GAG"]].upper()
                  and "X" not in seqs[domains["RT"]].upper())
    print("elements with GAG and RT, both X-free: %d" % len(keep))

    plaac = {r["SEQid"].strip(): r for r in
             csv.DictReader(open(os.path.join(GYDB, "plaac_data", "all_cores_a05.tsv"),
                                 encoding="utf-8"), delimiter="\t")}
    copia = {line.split("\t")[0] for line in
             open(os.path.join(GYDB, "copia_tree_tips.txt"), encoding="utf-8")
             if line.strip()}

    for kind, protein in (("rt_272.faa", "RT"), ("gag_272.faa", "GAG")):
        with open(os.path.join(GYDB, kind), "w", encoding="utf-8", newline="\n") as fh:
            for name in keep:
                fh.write(">%s\n%s\n" % (tip_label(name), seqs[elements[name][protein]]))

    table = read_superfamilies(SUPERFAMILY_CSV)
    lowered = {k.lower(): k for k in table}

    with open(os.path.join(GYDB, "traits_272.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["tip", "element", "group", "superfamily", "host", "clade",
                         "on_copia_tree", "rt_len", "short_rt", "gag_len", "llr",
                         "has_prd", "prd_start", "prd_end", "prd_mid_pct", "prd_zone",
                         "cchc", "cchc_n"])
        for name in keep:
            row = plaac[elements[name]["GAG"]]
            gag_seq = seqs[elements[name]["GAG"]]
            knuckles = CCHC.findall(gag_seq)
            gag_len, rt_len = int(row["PROTlen"]), len(seqs[elements[name]["RT"]])
            has_prd = int(row["PRDlen"]) > 0
            start, end = int(row["PRDstart"]), int(row["PRDend"])
            mid = 100.0 * (start + end) / 2 / gag_len if has_prd else None
            info = classify(name, table, lowered)
            writer.writerow([
                tip_label(name), name, group_of(name), info["superfamily"],
                info["host"], info["clade"],
                int(name in copia), rt_len, int(rt_len < SHORT_RT), gag_len,
                row["LLR"], int(has_prd),
                start if has_prd else "", end if has_prd else "",
                "%.1f" % mid if has_prd else "", prd_zone(mid),
                int(bool(knuckles)), len(knuckles)])

    families = collections.Counter(
        classify(n, table, lowered)["superfamily"] for n in keep)
    zones = collections.Counter(
        prd_zone(100.0 * (int(plaac[elements[n]["GAG"]]["PRDstart"]) +
                          int(plaac[elements[n]["GAG"]]["PRDend"])) / 2 /
                 int(plaac[elements[n]["GAG"]]["PROTlen"]))
        for n in keep if int(plaac[elements[n]["GAG"]]["PRDlen"]) > 0)
    print("  superfamilies: %s" % dict(families))
    print("  PrLD position: %s" % dict(zones))

    short = sum(1 for n in keep if len(seqs[elements[n]["RT"]]) < SHORT_RT)
    hits = sum(1 for n in keep if int(plaac[elements[n]["GAG"]]["PRDlen"]) > 0)
    groups = collections.Counter(group_of(n) for n in keep)
    print("  groups: %s" % dict(groups))
    print("  RT < %d aa (flagged, kept): %d" % (SHORT_RT, short))
    print("  PrLD called: %d" % hits)
    print("  on the Ty1/copia tree: %d" % len(set(keep) & copia))


if __name__ == "__main__":
    main()
