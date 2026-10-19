#!/usr/bin/env python3
"""Compare genodiff called VCF vs simulator truth VCF, with indel
normalization: indels are compared as (chrom, start, deleted_seq, inserted_seq)
after applying them to the reference, so any equivalent representation
(anchor base, right-expanded, etc.) compares equal.

Usage:
    python compare_vcf_norm.py truth.vcf called.vcf reference.fa
"""
import sys

def load(path):
    snps = set()
    indels = {}  # (chrom, start) -> (deleted, inserted)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, _, _, info = line.split("\t")[:8]
            pos = int(pos)
            if alt.startswith("<"):
                continue
            if len(ref) == 1 and len(alt) == 1:
                snps.add((chrom, pos, ref, alt))
            else:
                # VCF anchored form: POS is the anchor base (1-based),
                # REF/ALT share it
                start = pos  # 0-based position of the anchor base
                deleted = ref[1:]
                inserted = alt[1:]
                # trim common suffix, then common prefix
                while deleted and inserted and deleted[-1] == inserted[-1]:
                    deleted, inserted = deleted[:-1], inserted[:-1]
                while deleted and inserted and deleted[0] == inserted[0]:
                    deleted, inserted = deleted[1:], inserted[1:]
                    start += 1
                indels[(chrom, start)] = (deleted, inserted)
    return snps, indels

def load_ref(ref_fasta):
    refs = {}
    name = None
    buf = []
    with open(ref_fasta) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    refs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name:
            refs[name] = "".join(buf)
    return refs

def left_roll(indels, refs):
    """Roll every indel to its leftmost equivalent position against the
    reference, so any two equivalent representations compare equal."""
    out = {}
    for (chrom, start), (dele, ins) in indels.items():
        rseq = refs.get(chrom)
        if rseq is None:
            out[(chrom, start)] = (dele, ins)
            continue
        s = start
        # deletion: can shift left while last base of deleted == base before
        while s > 0 and dele and rseq[s + len(dele) - 1] == rseq[s - 1]:
            dele = rseq[s - 1] + dele[:-1]
            s -= 1
        # insertion: can shift left while last base of inserted == base before
        while s > 0 and ins and rseq[s - 1] == ins[-1]:
            ins = ins[-1] + ins[:-1]
            s -= 1
        out[(chrom, s)] = (dele, ins)
    return out

def main(truth_vcf, called_vcf, ref_fasta):
    refs = load_ref(ref_fasta)
    t_snps, t_ind = load(truth_vcf)
    c_snps, c_ind = load(called_vcf)
    t_ind = left_roll(t_ind, refs)
    c_ind = left_roll(c_ind, refs)

    tp = len(t_snps & c_snps)
    print(f"SNP    : truth={len(t_snps)} called={len(c_snps)} TP={tp} "
          f"FP={len(c_snps - t_snps)} FN={len(t_snps - c_snps)} "
          f"(P={tp/max(1,len(c_snps))*100:.2f}% R={tp/max(1,len(t_snps))*100:.2f}%)")

    matched = 0
    fp, fn = [], []
    for key, val in c_ind.items():
        if t_ind.get(key) == val:
            matched += 1
        else:
            fp.append((key, val, t_ind.get(key)))
    for key, val in t_ind.items():
        if c_ind.get(key) != val:
            fn.append((key, val, c_ind.get(key)))
    print(f"InDel  : truth={len(t_ind)} called={len(c_ind)} TP={matched} "
          f"FP={len(fp)} FN={len(fn)} "
          f"(P={matched/max(1,len(c_ind))*100:.2f}% R={matched/max(1,len(t_ind))*100:.2f}%)")
    for tag, items in [("FP", fp[:5]), ("FN", fn[:5])]:
        for it in items:
            print(f"  {tag}: at {it[0]} called={str(it[1])[:50]} truth={str(it[2])[:50]}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
