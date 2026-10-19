#!/usr/bin/env python3
"""Union-merge single-sample VCFs into ONE multi-sample VCF.

Semantics: every (CHROM, POS, REF, ALT, INFO) record seen in ANY input becomes
one row; a sample's GT is 1 if it carries that exact record, 0 (reference) if
it does not. Haploid calling convention. SV annotations in INFO (SVTYPE, SVLEN,
CHR2, POS2) are carried through, so structural variants stay self-describing
after the merge.

Rows are sorted by contig (in header order), then position, then alleles.
Usage:
    python merge_vcfs.py -o all_samples.vcf s1.vcf s2.vcf [...]
"""
import argparse
import os
import sys

# declared once in the merged header so symbolic alleles and the SV INFO
# fields are self-describing
SV_HEADERS = [
    '##ALT=<ID=DEL,Description="Deletion">',
    '##ALT=<ID=INS,Description="Insertion">',
    '##ALT=<ID=INV,Description="Inversion">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=BND,Description="Breakend">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the structural variant">',
    '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome of second breakend">',
    '##INFO=<ID=POS2,Number=1,Type=Integer,Description="Position of second breakend">',
    '##INFO=<ID=QSEQ,Number=1,Type=String,Description="Query contig of this variant\'s query-side span (per-sample files in ##genodiff_query lines)">',
    '##INFO=<ID=QSTART,Number=1,Type=Integer,Description="1-based start of the query-side span; extract with e.g. samtools faidx query.fa QSEQ:QSTART-QEND">',
    '##INFO=<ID=QEND,Number=1,Type=Integer,Description="1-based end of the query-side span">',
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcfs", nargs="+", help="one or more single-sample VCFs")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    # exclude the output file and any multi-sample VCFs from the inputs
    out_abs = os.path.abspath(a.out)
    inputs = []
    for p in a.vcfs:
        if os.path.abspath(p) == out_abs:
            continue
        # skip files that are already multi-sample (contain >1 sample column)
        with open(p) as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    n_samples = len(line.rstrip().split("\t")) - 9
                    if n_samples > 1:
                        print(f"  SKIP {p}: already multi-sample ({n_samples} samples)",
                              file=sys.stderr)
                    else:
                        inputs.append(p)
                    break
    if not inputs:
        sys.exit("no single-sample VCFs to merge")
    a.vcfs = inputs

    samples = []
    contigs = []  # (name, full ##contig line); first definition wins, unless a
                  # later one adds the length field
    reference = ""      # ##reference of the first input (all inputs share it)
    query_files = []    # ##genodiff_query lines, one per input sample
    masks = {}  # "chrom\tpos\tref\talt\tinfo" -> bitmask over samples

    for si, path in enumerate(a.vcfs):
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith("##contig=<ID="):
                    name = line.split("ID=")[1].split(",")[0].rstrip(">")
                    for i, (n, cur) in enumerate(contigs):
                        if n == name:
                            if "length=" in line and "length=" not in cur:
                                contigs[i] = (name, line)
                            break
                    else:
                        contigs.append((name, line))
                    continue
                if line.startswith("##reference=") and not reference:
                    reference = line
                    continue
                if line.startswith("##genodiff_query="):
                    if line not in query_files:
                        query_files.append(line)
                    continue
                if line.startswith("#CHROM"):
                    samples.append(line.split("\t")[-1])
                    continue
                if line.startswith("#"):
                    continue
                f = line.split("\t", 8)
                if len(f) < 5:
                    continue
                chrom, pos, ref, alt, info = f[0], f[1], f[3], f[4], f[7]
                key = f"{chrom}\t{pos}\t{ref}\t{alt}\t{info}"
                bit = 1 << si
                if key in masks:
                    masks[key] |= bit
                else:
                    masks[key] = bit
        print(f"  read {path}: sample={samples[-1]}", file=sys.stderr)

    if len(set(samples)) != len(samples):
        sys.exit("duplicate sample names in input VCFs")

    rank = {n: i for i, (n, _) in enumerate(contigs)}

    def sort_key(k):
        chrom, pos, ref, alt, info = k.split("\t")
        return (rank.get(chrom, 1 << 30), chrom, int(pos), ref, alt, info)

    keys = sorted(masks.keys(), key=sort_key)

    with open(a.out, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=genodiff-merge\n")
        if reference:
            out.write(reference + "\n")
        for q in query_files:
            out.write(q + "\n")
        for _, c in contigs:
            out.write(c + "\n")
        for h in SV_HEADERS:
            out.write(h + "\n")
        out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype (0=ref, 1=alt; haploid)\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples) + "\n")

        for key in keys:
            chrom, pos, ref, alt, info = key.split("\t")
            bits = masks[key]
            gts = "\t".join("1" if (bits >> si) & 1 else "0" for si in range(len(samples)))
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\tGT\t{gts}\n")

    print(f"wrote {a.out}: {len(keys)} sites x {len(samples)} samples", file=sys.stderr)


if __name__ == "__main__":
    main()
