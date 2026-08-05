import csv
import os
from collections import defaultdict

REXDB = os.path.join("data", "REXdb")
OUT = "data"
PROTEIN_FILE = "Viridiplantae_v4.0.fasta"
CLASSIFICATION_FILE = "Viridiplantae_v4.0.classification"
INFO_FILE = "Viridiplantae_v4.0_ALL_info-species-source_taxid"
DNA_FILE = "Viridiplantae_v4.0_ALL_DNA.fasta"

SUPERFAMILIES = {"Ty1/copia": "copia", "Ty3/gypsy": "gypsy"}
UPSTREAM = 75
GAG_MIN_UPSTREAM = 1

# S. cerevisiae Ty1/copia elements, also kept at full length as controls (F-10).
REFERENCES = {"REXdb_ID5405", "REXdb_ID5406", "REXdb_ID5412",
              "REXdb_ID5413", "REXdb_ID5415"}

# The standard genetic code. "*" is a stop codon -- that is what the Gag
# extension is cut at, and what tells us a reading frame has ended.
CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",

    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",

    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# gag_status says why a trait is missing, so absence is never ambiguous (F-7):
#   ok               extracted and written to gag_*.faa
#   no_slice         REXdb published no Gag for this element -- genuine absence
#   unplaceable      slice exists but is not contiguous in its own DNA, usually
#                    a frameshift -- the element keeps its place in the tree
#   record_mismatch  no part of the slice is in the DNA at all, so the two
#                    release files disagree; dropped from the tree too (F-16)
#   no_nterm         placed, but nothing upstream of the core (F-8)
COLUMNS = ["rexdb_id", "superfamily", "lineage", "species", "taxid", "prelim_id",
           "source", "is_reference", "rt_len", "rt_ambiguous", "in_tree",
           "has_gag", "gag_status", "gag_core_exact", "gag_frame",
           "gag_upstream", "gag_core_len", "gag_len", "gag_ambiguous",
           "gag_clipped_by"]


def element_number(element):
    """REXdb_ID4903 -> 4903, so output sorts numerically."""
    return int(element.replace("REXdb_ID", ""))


def translate_codons(dna):
    """Read three bases at a time into amino acids. Unexpected bases give X."""
    return "".join(CODONS.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))


def forward_frames(dna):
    """Yield (name, protein) for the three forward reading frames.

    Reverse frames are not searched: REXdb DNA records are already oriented, and
    all 12,806 located elements matched on a forward frame (F-13).
    """
    for offset in (0, 1, 2):
        yield "+%d" % (offset + 1), translate_codons(dna[offset:])


def search_frames(dna, query):
    """Find query in a forward frame, X on either side matching anything.

    REXdb writes X where it could not resolve a codon, and str.find cannot match
    a wildcard. So the query is anchored on its longest X-free run and only those
    positions are compared -- checking every position would be far too slow. A
    query with no X anchors on itself, which makes this an exact search.
    """
    anchor = max(query.split("X"), key=len)
    if not anchor:
        return None, None, -1
    offset = query.find(anchor)
    for frame, protein in forward_frames(dna):
        start = 0
        while (index := protein.find(anchor, start)) >= 0:
            begin = index - offset
            window = protein[begin:begin + len(query)]
            if begin >= 0 and all(a == b or "X" in (a, b)
                                  for a, b in zip(window, query)):
                return frame, protein, begin
            start = index + 1
    return None, None, -1


def find_slice(dna, slice_):
    """Locate a published slice. Returns (frame, protein, index, exact)."""
    frame, protein, index = search_frames(dna, slice_)
    if frame:
        return frame, protein, index, True
    for head in (30, 25, 20):
        if len(slice_) <= head:
            break
        frame, protein, index = search_frames(dna, slice_[:head])
        if frame:
            return frame, protein, index, False
    return None, None, -1, False


def slice_absent(dna, slice_):
    """True when no 20 aa stretch of the slice appears in any frame.

    A frameshifted Gag still has most of itself in the DNA, just not contiguously.
    Nothing matching at all means the protein and DNA records describe different
    sequences, so the element is untrustworthy and leaves the tree too (F-16).
    """
    frames = [protein for _, protein in forward_frames(dna)]
    return not any(slice_[i:i + 20] in protein
                   for i in range(0, len(slice_) - 19, 10)
                   for protein in frames)


def extend_upstream(protein, start, limit):
    """Sequence before `start`, at most `limit` aa (None = all), cut at any
    in-frame stop. Returns (sequence, reason the extension ended)."""
    begin = 0 if limit is None else max(0, start - limit)
    region = protein[begin:start]
    stop = region.rfind("*")
    if stop >= 0:
        return region[stop + 1:], "stop"
    return region, "" if begin > 0 else "record_start"


def read_fasta(path):
    """Yield (header, sequence)."""
    header, parts = None, []
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.strip()
        if line.startswith(">"):
            if header:
                yield header, "".join(parts)
            header, parts = line[1:].strip(), []
        elif line:
            parts.append(line)
    if header:
        yield header, "".join(parts)


def read_metadata():
    """{rexdb_id: row} for copia and gypsy elements only."""
    classification = {}
    with open(os.path.join(REXDB, CLASSIFICATION_FILE), encoding="utf-8",
              errors="replace") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            fields = [f.strip() for f in fields]
            if fields and fields[0]:
                classification[fields[0]] = (fields[3] if len(fields) > 3 else "",
                                             fields[4] if len(fields) > 4 else "")

    metadata = {}
    with open(os.path.join(REXDB, INFO_FILE), encoding="utf-8",
              errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            superfamily, lineage = classification.get(row["REXdb_ID"], ("", ""))
            if superfamily in SUPERFAMILIES:
                metadata[row["REXdb_ID"]] = dict(
                    rexdb_id=row["REXdb_ID"], superfamily=superfamily,
                    lineage=lineage,
                    species=" ".join(row["Species"].split()),  # stray double spaces
                    taxid=row["NCBI_TaxID"], prelim_id=row["prelim.ID"],
                    source=row["Data.source"],
                    is_reference=int(row["REXdb_ID"] in REFERENCES))
    return metadata


def read_domains(wanted):
    """{rexdb_id: {domain: sequence}} for RT and GAG only."""
    domains = defaultdict(dict)
    for header, sequence in read_fasta(os.path.join(REXDB, PROTEIN_FILE)):
        left, _, element = header.partition("__")
        domain = left.rsplit("-", 1)[-1]
        if element in wanted and domain in ("RT", "GAG"):
            domains[element][domain] = sequence.upper()
    return domains


def write_fasta(name, records):
    with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="\n") as handle:
        for element, sequence in sorted(records, key=lambda r: element_number(r[0])):
            handle.write(">%s\n%s\n" % (element, sequence))


def main():
    print("running...")
    rows = read_metadata()
    domains = read_domains(set(rows))

    for element, row in rows.items():
        rt = domains[element].get("RT", "")
        has_gag = "GAG" in domains[element]
        row.update(rt_len=len(rt), rt_ambiguous=rt.count("X"),
                   in_tree=int(bool(rt)), has_gag=int(has_gag),
                   gag_core_len=len(domains[element].get("GAG", "")),
                   gag_status="unplaceable" if has_gag else "no_slice",
                   gag_core_exact="", gag_frame="", gag_upstream="",
                   gag_len="", gag_ambiguous="", gag_clipped_by="")

    gag = defaultdict(list)
    references = []
    for element, dna in read_fasta(os.path.join(REXDB, DNA_FILE)):
        core = domains.get(element, {}).get("GAG")
        if not core:
            continue
        sequence_dna = dna.upper()
        frame, protein, index, exact = find_slice(sequence_dna, core)
        if not frame:
            if slice_absent(sequence_dna, core):
                rows[element].update(gag_status="record_mismatch", in_tree=0)
            continue
        upstream, ended_by = extend_upstream(protein, index, UPSTREAM)
        sequence = upstream + core
        # Gag needs a usable RT too: an element absent from the tree has nowhere
        # to carry a trait.
        usable = bool(rows[element]["rt_len"]) and len(upstream) >= GAG_MIN_UPSTREAM
        rows[element].update(gag_status="ok" if usable else "no_nterm",
                             gag_core_exact=int(exact), gag_frame=frame,
                             gag_upstream=len(upstream), gag_len=len(sequence),
                             gag_ambiguous=sequence.count("X"),
                             gag_clipped_by=ended_by)
        if usable:
            gag[rows[element]["superfamily"]].append((element, sequence))
        if element in REFERENCES:  # controls bypass the filters
            full, _ = extend_upstream(protein, index, None)
            references.append((element, full + core))

    for superfamily, short in SUPERFAMILIES.items():
        write_fasta("rt_%s.faa" % short,
                    [(e, domains[e]["RT"]) for e, row in rows.items()
                     if row["superfamily"] == superfamily and row["in_tree"]])
        write_fasta("gag_%s.faa" % short, gag[superfamily])
    write_fasta("gag_reference.faa", references)

    with open(os.path.join(OUT, "elements.tsv"), "w", encoding="utf-8",
              newline="\n") as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for element in sorted(rows, key=element_number):
            writer.writerow(rows[element])


if __name__ == "__main__":
    main()
