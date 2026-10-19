# genodiff

[![Rust](https://img.shields.io/badge/language-Rust-DEA584.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-45%20passed-brightgreen.svg)](#testing)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20Windows-lightgrey.svg)](#installation)

**A self-contained, zero-dependency tool for whole-genome alignment and variant
calling between haploid genome assemblies.**

genodiff takes two (or more) assembled genome FASTA files and produces a
standard VCF 4.2 file containing SNPs, small insertions/deletions, and
structural variants (DEL / INS / INV / DUP / BND) — all from a single static
Rust binary with **no external dependencies**.

> 📖 For the full pipeline description, algorithm details, and innovation
> points, see [docs/PIPELINE.md](docs/PIPELINE.md).
> 📖 中文文档见 [README_CN.md](README_CN.md)。

---

## Why genodiff?

| | genodiff | minimap2 + Syri | MUM&Co |
|---|---|---|---|
| External dependencies | **0** | ≥4 (minimap2, Python3, NumPy, SciPy, Matplotlib, NetworkX, Pysam) | MUMmer 3/4 + Python |
| One command FASTA→VCF | **Yes** (`call`) | Multi-step pipeline | Partial |
| SV types | DEL/INS/INV/DUP/BND | DEL/INS/INV/DUP/TRA | Partial |
| Direct VCF output (left-aligned) | **Yes** | Custom + VCF | VCF |
| Built-in variant simulator + truth VCF | **Yes** | No | No |
| Built-in multi-sample merge | **Yes** | Needs bcftools | No |
| Memory model | Sorted array, ~24 B/minimizer | Index + SAM intermediate | Delta intermediate |
| Install | `cargo build` (no deps to download) | Multi-tool install | MUMmer build |

### Key innovations

1. **Single static binary, zero dependencies** — one Rust binary goes from
   FASTA input to VCF output with no interpreter, no third-party crates,
   and no intermediate SAM/PAF conversion layer.
2. **Purpose-built for assembly-to-assembly comparison** — bypasses the
   read-mapping parameter space entirely; collinear anchor chains are the
   primary data structure (chain = alignment block = SV candidate).
3. **Patience-based chain tolerance** — keeps chains open through
   mid-frequency repeats that would otherwise shatter them into 4–5-anchor
   fragments (on real 800 Mb plant genomes: 1.03×10⁷ spurious breaks →
   5.2×10⁵, anchors in chains +50%).
4. **Symmetric-gap inversion calling** — inversion flanks separated by
   shared repeats of any length are bridged when query-side and reference-side
   gap lengths match, overcoming fixed-tolerance limits.
5. **Seed-overlap slip indel recovery** — detects indels hidden inside
   overlapping seed claims (minimizers are 19 bp but sampled ~every 10 bp),
   recovering homopolymer and short-repeat indels missed by conventional
   gap alignment.
6. **Repeat-attractor filtering** — automatically identifies and excludes
   paralogous pileups (16 kb reference bins hit by ≥5 mutually distant query
   regions), suppressing false SNPs and SVs in repeat-dense regions.
7. **Bounded-memory streaming design** — index is a single sorted array
   (~24 B/minimizer), chain state is O(1), contig-level parallelism with
   in-flight byte budget; 800 Mb × 800 Mb full pipeline in 117 s / ~4 GB
   (48 threads).
8. **Built-in validation loop** — deterministic variant simulator producing
   truth VCFs + left-roll normalization comparison script; on simulated
   genomes: SNP/InDel precision & recall 100%, all four SV classes 100%
   recalled; 45 unit/integration tests covering every algorithm stage.

---

## Installation

```bash
# Install Rust (skip if already installed; user-level, no root needed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"

git clone https://github.com/<your-username>/genodiff.git
cd genodiff
cargo build --release          # binary at target/release/genodiff

# Verify: 45 tests should pass
cargo test --release
```

No third-party crates are used; compilation requires only the Rust toolchain.

## Quick start

```bash
G=target/release/genodiff

# One command: FASTA in → VCF out (SNP + InDel + SV)
$G call reference.fa query.fa -o query_vs_ref.vcf --threads 48 --sample query_name
```

That's it. Internally, `call` builds the reference index in memory, chains
seed anchors, banded-aligns the gaps, calls SNPs/InDels, infers SVs from
chain topology, and writes a standard VCF 4.2 file.

## Step-by-step (PAF + reusable index)

```bash
# ① Build a reference index once (.gdx file, reusable for all queries)
$G index reference.fa reference.gdx         # k=19, w=19 by default

# ② Chain-align any number of queries → PAF alignment blocks
$G map reference.gdx query.fa --threads 48 2>map.log > query_vs_ref.paf
```

PAF output (one line per collinear block):
`q_name q_len q_start q_end strand r_name r_len r_start r_end n_anchors na:n_anchors`

## Multi-genome batch mode

Edit the `QUERIES` array in `scripts/run_all.sh` (label=path pairs), then:

```bash
bash scripts/run_all.sh
# Produces: batch/<sample>.vcf (one per genome) + batch/all_samples.vcf (merged)
```

The merge produces a multi-sample VCF where every variant site seen in any
sample gets one row; GT=1 if the sample carries it, GT=0 if reference.

## Visualization (PAF → publication-quality figures)

Three pure-matplotlib scripts (requires Python3 + numpy + matplotlib):

```bash
# Dot plot: reference along x, query along y
python scripts/plot_syntenic_dotplot.py query_vs_ref.paf -o dotplot.png --top 8

# Linear synteny: reference row (top), query row (bottom), ribbon-connected
python scripts/plot_syntenic_linear.py query_vs_ref.paf -o synteny.svg --top 8

# Circos: reference semicircle (top), query semicircle (bottom)
python scripts/plot_syntenic_circos.py query_vs_ref.paf -o circos.png --top 8
```

All support `.png`, `.pdf`, and hybrid `.svg` output (900-dpi raster base +
individually editable `<text>` elements, pixel-exact alignment).

## Validation

```bash
# Generate a known-mutant genome + truth VCF
$G sim reference.fa simdir --snp 500 --indel 50 --svdel 5 --svins 5 --inv 2 --dup 2

# Call variants and compare
$G call reference.fa simdir/query.fa -o called.vcf
python scripts/compare_vcf_norm.py simdir/truth.vcf called.vcf reference.fa
# Expected: SNP P=100% R=100%; InDel P=100% R=100%
```

On a 5 Mb genome with 5,530 injected variants (5,000 SNP + 500 small InDel +
10 DEL + 10 INS + 5 INV + 5 DUP): **SNP 100%/100%, InDel 100%/100%,
all four SV classes 100% recalled**.

On real data (maize zm1 × caerulea, ~800 Mb each, 48 threads):
**117 s full pipeline** producing 5.97M SNPs, 1.24M+1.74M InDels, and
55K SVs (DEL=25,228 INS=17,408 INV=10 DUP=1,831 BND=5,056).

## Algorithm overview

```
FASTA input
  → ② minimizer index (k=19, w=19, canonical, sorted array)
  → ③ seed anchors (binary search + max-occ repeat filter)
  → ④ greedy chaining (drift tolerance + patience + short/secondary filters)
       → PAF alignment blocks
  → ⑤ chain-topology SV inference (DEL/INS/INV/DUP/BND + repeat filter)
  → ⑥ banded affine alignment (Gotoh: match +2, mismatch −6, open −8, ext −2)
  → ⑦ variant extraction (SNP/InDel + left-normalization + strand correction)
  → ⑧ VCF 4.2 output → multi-sample merge + visualization
```

See [docs/PIPELINE.md](docs/PIPELINE.md) for the full 8-stage description
with per-stage parameters, innovation details, and a Mermaid flowchart.

## Input compatibility

- Haploid assemblies (any number of contigs); ref/query roles are swappable
- FASTA headers: first whitespace-delimited token is used as sequence name
  (handles trailing spaces, descriptions)
- Upper/lowercase bases unified; N and IUPAC ambiguity codes skipped
- Windows/Unix line endings; empty lines tolerated
- Contig names need not correspond between the two genomes

## CLI reference

```
genodiff index   <ref.fa> <out.gdx> [--k 19] [--w 19]
genodiff map     <ref.gdx> <query.fa> [--threads N] [--max-occ 100]
                 [--drift 150] [--max-gap 100000] [--min-anchors 25]
                 [--min-span 500] [--patience 1000] [--keep-secondary]
genodiff call    <ref.fa> <query.fa> [-o out.vcf] [--sample NAME]
                 [--threads N] [--max-occ 100] [--drift 150]
                 [--max-gap 100000] [--min-anchors 25] [--min-span 500]
                 [--patience 1000] [--keep-secondary]
genodiff sim     <ref.fa> <outdir> [--seed 42] [--snp 500] [--indel 50]
                 [--svdel 5] [--svins 5] [--inv 2] [--dup 2]
genodiff randref <out.fa> <length> [--contigs 1] [--seed 1]
```

## Known limitations

1. SV breakpoint precision is chain-end level (~hundreds of bp); local
   realignment for base-level breakpoints is planned.
2. Greedy chaining may select a suboptimal paralog copy in multi-copy
   regions (DP chaining is on the roadmap).
3. Complex rearrangements with both sides >2 kb apart and no shared repeats
   are reported as chain breaks without further classification.
4. Centromeric / highly repetitive regions are unalignable at the seed
   level (a fundamental limitation shared by all sequence-based methods).

## Roadmap

- Resident index server (eliminate redundant batch re-indexing)
- DP chaining (O(n log n)) for fragmented assemblies
- SV breakpoint refinement (local realignment to base-pair precision)
- Multi-sample VCF phasing and consensus generation

## License

[MIT](LICENSE)
