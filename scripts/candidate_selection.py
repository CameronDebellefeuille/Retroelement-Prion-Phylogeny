"""Pull the shortlisted Gag candidates out of GyDB and write them up for synthesis.

Fourteen elements in three groups. Eleven are PrLD-bearing: six Ty3/Gypsy from
insects, three Retroviridae foamy viruses, and two Ty1/Copia yeast elements.
The last three are high-disorder controls with no PrLD at all. Edit CANDIDATES
to change the list -- which group an element lands in is read from has_prd in
traits_272.tsv, not hardcoded here, so the two cannot drift apart.

The no-PrLD three are the disorder control: they answer whether the PrLD is
doing something an ordinary IDR of the same Gag does not. They are the three
highest mean_disorder elements with has_prd == 0 and a Gag of at least 60 aa.
The 60 aa floor is PLAAC's core window -- below it PLAAC cannot call a domain
at all, so "no PrLD" would be absence of evidence rather than evidence of
absence. Two CopiaSL_monotypic fragments of 26 and 25 aa score higher on raw
disorder and are excluded for exactly that reason.

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
INTERPRO = os.path.join(REPO, "data", "processed", "gydb", "InterProScan")
DISORDER = os.path.join(REPO, "data", "processed", "gydb", "disorder_272.tsv")
OUT = os.path.join(REPO, "candidates")
LLPS = os.path.join(OUT, "scores.csv")

# GyDB has no element called "Ty1". The S. cerevisiae Ty1 element is
# catalogued as Ty1B, which is what goes in the list and what comes back out in
# the element column -- renaming it here would break the join to traits_272.tsv
# and to GAG_Ty1B in the cores database.
CANDIDATES = ["Tom", "297", "17.6", "Yoyo", "HMS-Beagle", "Nomad",
              "EFV", "BFV", "FFV",
              "Ty1B", "Tkm1",
              "CopiaSL_29", "RIRE2", "RetroSor1"]

# A zinc knuckle is a hit to any of these four signatures for IPR001878 /
# IPR036875 -- Pfam, SMART, PROSITE and SUPERFAMILY. Same definition the knuckle
# figures use, so the calls here and there cannot drift apart.
#
# NOT the cchc column in traits_272.tsv. That is a regex, and nc_scan.py
# measured it calling a knuckle in 70% of copia Gags where the profile models
# call 21% -- a 14-residue pattern with four fixed positions hits by chance in a
# 300 aa protein.
KNUCKLE = ("PF00098", "SM00343", "PS50158", "SSF57756")


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


def composition(seq):
    """Residue percentages, net charge and charged fraction.

    Q and N are reported together as well as apart: Q/N richness is the classic
    prion signature and what PLAAC's HMM is largely responding to, so the pair
    is the number to compare across candidates. G, S, Y and P are the other
    residues that recur in characterised PrLDs -- Y in particular, because pi-pi
    stacking is a proposed driver of the phase behaviour.

    Net charge is at neutral pH, counting K and R as +1 and D and E as -1;
    histidine is left out rather than counted as a fraction. Low net charge with
    high Q/N is the composition that separates a prion-like domain from an
    ordinary disordered region, which tends to be charge-driven.
    """
    n = len(seq)
    if not n:
        return {}
    pct = lambda residues: round(100 * sum(seq.count(r) for r in residues) / n, 1)
    pos, neg = sum(seq.count(r) for r in "KR"), sum(seq.count(r) for r in "DE")
    return {"pct_Q": pct("Q"), "pct_N": pct("N"), "pct_QN": pct("QN"),
            "pct_G": pct("G"), "pct_S": pct("S"), "pct_Y": pct("Y"),
            "pct_P": pct("P"), "net_charge": pos - neg,
            "pct_charged": round(100 * (pos + neg) / n, 1)}


def read_disorder(path):
    """metapredict over the whole Gag core, from disorder_272.tsv.

    Whole Gag, NOT the PrLD on its own -- disorder_272.py scores each core in
    one pass and the table keeps no per-region breakdown. Recomputing over the
    domain alone would mean rerunning metapredict, which is the one dependency
    not installed by default. So read these as context for the protein the
    domain sits in, not as a measurement of the domain.
    """
    return {r["tip"]: r for r in
            csv.DictReader(open(path, encoding="utf-8"), delimiter="\t")}


def read_llps(path):
    """catGRANULE 2.0 LLPS propensity, one score per candidate.

    Run externally on the web service -- L-7 leaves catGRANULE as an external
    input, so scores.csv is pasted in rather than computed here. It was run on
    candidate_prld.faa, i.e. the DOMAIN, not the whole Gag: every profile in
    data.csv is exactly prld_len long. So the score is the propensity of the
    prion-like domain on its own, which is not the same question as whether the
    intact Gag phase separates.

    Missing file is not an error -- the table is still useful without it.
    """
    if not os.path.exists(path):
        return {}
    out = {}
    for r in csv.DictReader(open(path, encoding="utf-8-sig")):
        name = (r.get("Name") or "").strip()
        if name:
            out[name] = r["LLPS_Score"].strip()
    return out


def read_interpro(directory):
    """Return (scanned elements, knuckle hits per element).

    Both are needed, and conflating them is the mistake to avoid: an element
    with no knuckle row and an element that was never scanned look identical if
    you only ask "any hits?". The first is evidence of absence, the second is
    absence of evidence.
    """
    scanned, hits = set(), {}
    for name in sorted(os.listdir(directory)):
        if not name.startswith("part"):
            continue
        for line in open(os.path.join(directory, name), encoding="utf-8",
                         errors="replace"):
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            scanned.add(f[0].strip())
            if f[4] in KNUCKLE:
                hits.setdefault(f[0].strip(), []).append((f[4], f[6], f[7]))
    return scanned, hits


def main():
    seqs = read_cores(CORES)
    traits = {r["element"]: r for r in
              csv.DictReader(open(TRAITS, encoding="utf-8"), delimiter="\t")}

    missing = [c for c in CANDIDATES if c not in traits]
    if missing:
        sys.exit("not in traits_272.tsv: %s" % ", ".join(missing))

    scanned, knuckles = read_interpro(INTERPRO)
    disorder = read_disorder(DISORDER)
    llps = read_llps(LLPS)

    rows = []
    for name in CANDIDATES:
        t = traits[name]
        gag = seqs["GAG_" + name]
        # The no-PrLD controls have empty prd_start/prd_end, so there is no
        # domain to slice and every prld_* column stays blank. Blank, not zero:
        # a zero would read as a measured composition of 0%, which is a
        # different claim from "there is no domain to measure".
        has_prld = t["has_prd"] == "1"
        if has_prld:
            # PLAAC coordinates are 1-based and inclusive.
            start, end = int(t["prd_start"]), int(t["prd_end"])
            prld = gag[start - 1:end]
        else:
            start, end, prld = "", "", ""
        hit = knuckles.get(name, [])
        if name not in scanned:
            knuckle = "not scanned"
        else:
            knuckle = "yes" if hit else "no"
        pc = composition(prld)
        gc = composition(gag)
        d = disorder.get(t["tip"], {})
        rows.append({
            "element": name,
            "superfamily": t["superfamily"],
            "host": t["host"],
            "gag_len": len(gag),
            "has_prld": int(has_prld),
            "prld_start": start,
            "prld_end": end,
            "prld_len": len(prld) if has_prld else "",
            "plaac_llr": t["llr"],
            "zinc_knuckle": knuckle,
            "knuckle_n": len(hit),
            "prld_pct_Q": pc.get("pct_Q", ""),
            "prld_pct_N": pc.get("pct_N", ""),
            "prld_pct_QN": pc.get("pct_QN", ""),
            "prld_pct_G": pc.get("pct_G", ""),
            "prld_pct_S": pc.get("pct_S", ""),
            "prld_pct_Y": pc.get("pct_Y", ""),
            "prld_pct_P": pc.get("pct_P", ""),
            "prld_net_charge": pc.get("net_charge", ""),
            "prld_pct_charged": pc.get("pct_charged", ""),
            "gag_pct_QN": gc["pct_QN"],
            "gag_mean_disorder": d.get("mean_disorder", ""),
            "llps_score": llps.get(name, ""),
            "gag_aa": gag,
            "prld_aa": prld})

    os.makedirs(OUT, exist_ok=True)

    def fasta(path, key, label, records):
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("; %s -- PROTEIN (amino acid), not DNA. See "
                     "scripts/candidate_selection.py\n" % label)
            for r in records:
                if r["has_prld"]:
                    domain = ("PrLD %s-%s | LLR %s | zinc knuckle: %s | "
                              "Q+N %.1f%% | net charge %+d"
                              % (r["prld_start"], r["prld_end"], r["plaac_llr"],
                                 r["zinc_knuckle"], r["prld_pct_QN"],
                                 r["prld_net_charge"]))
                else:
                    domain = ("no PrLD (disorder control) | LLR %s | "
                              "zinc knuckle: %s | mean disorder %s"
                              % (r["plaac_llr"], r["zinc_knuckle"],
                                 r["gag_mean_disorder"] or "-"))
                fh.write(">%s | %s | %s | Gag %d aa | %s\n%s\n"
                         % (r["element"], r["superfamily"],
                            r["host"] or "host unrecorded", r["gag_len"],
                            domain, wrap(r[key])))

    # Every candidate has a Gag; only the eleven PrLD-bearing ones have a domain
    # to write, so candidate_prld.faa is the shorter file by design.
    fasta(os.path.join(OUT, "candidate_gag.faa"), "gag_aa", "full Gag cores", rows)
    fasta(os.path.join(OUT, "candidate_prld.faa"), "prld_aa",
          "PLAAC-called domains only", [r for r in rows if r["has_prld"]])

    cols = ["element", "superfamily", "host", "gag_len", "has_prld",
            "prld_start", "prld_end",
            "prld_len", "plaac_llr", "zinc_knuckle", "knuckle_n",
            "prld_pct_Q", "prld_pct_N", "prld_pct_QN", "prld_pct_G", "prld_pct_S",
            "prld_pct_Y", "prld_pct_P", "prld_net_charge", "prld_pct_charged",
            "gag_pct_QN", "gag_mean_disorder", "llps_score", "gag_aa", "prld_aa"]
    with open(os.path.join(OUT, "candidates.tsv"), "w", encoding="utf-8",
              newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    if llps:
        # Rank only the rows that have a score. The no-PrLD controls were never
        # sent to catGRANULE -- it was run on candidate_prld.faa, and they have
        # no domain in that file -- so they are absent here rather than ranked
        # last, which would read as a measured propensity of zero.
        scored = [r for r in rows if r["llps_score"]]
        ranked = sorted(scored, key=lambda r: -float(r["llps_score"]))
        lcols = ["rank", "element", "superfamily", "host", "llps_score",
                 "prld_len", "plaac_llr", "prld_pct_QN", "prld_net_charge",
                 "prld_pct_Y", "gag_mean_disorder", "zinc_knuckle"]
        with open(os.path.join(OUT, "LLPS_score.tsv"), "w", encoding="utf-8",
                  newline="\n") as fh:
            w = csv.DictWriter(fh, fieldnames=lcols, delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            for i, r in enumerate(ranked, 1):
                w.writerow({**{k: r[k] for k in lcols if k in r}, "rank": i})
        missing = [r["element"] for r in rows
                   if r["has_prld"] and not r["llps_score"]]
        if missing:
            print("no LLPS score for: %s" % ", ".join(missing))
    print("%d candidates -> %s" % (len(rows), OUT))
    print("%-12s %-13s %5s %10s %7s %8s %6s %6s %6s %9s" %
          ("element", "superfamily", "Gag", "PrLD", "LLR", "knuckle",
           "Q+N%", "charge", "Y%", "disorder"))
    for r in rows:
        if r["has_prld"]:
            dom = "%4d-%-5d" % (r["prld_start"], r["prld_end"])
            comp = ("%5.1f%% %+6d %5.1f%%"
                    % (r["prld_pct_QN"], r["prld_net_charge"], r["prld_pct_Y"]))
        else:
            dom, comp = "%10s" % "no PrLD", "%6s %6s %6s" % ("-", "-", "-")
        print("%-12s %-13s %5d %s %7.1f %8s %s %9s"
              % (r["element"], r["superfamily"], r["gag_len"], dom,
                 float(r["plaac_llr"]), r["zinc_knuckle"], comp,
                 r["gag_mean_disorder"] or "-"))

    called = sum(1 for r in rows if r["zinc_knuckle"] == "yes")
    unscanned = sum(1 for r in rows if r["zinc_knuckle"] == "not scanned")
    print("\n%d of %d carry a zinc knuckle (%d not scanned by InterProScan)."
          % (called, len(rows), unscanned))
    print("Knuckle = IPR001878 / IPR036875 via %s." % ", ".join(KNUCKLE))
    print("\nProtein sequences. For DNA, codon-optimise or pull the native "
          "nucleotide from GenBank -- see the module docstring.")


if __name__ == "__main__":
    main()
