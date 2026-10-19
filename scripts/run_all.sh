#!/bin/bash
# Batch mode: align every genome to zm1 reference ONE AT A TIME (sequential:
# each VCF is finished before the next genome starts, so peak disk usage is
# just one VCF at a time), then merge everything into a single multi-sample
# VCF at the end.
#
#   bash run_all.sh
#
# Output: batch/<name>.vcf (per genome) + batch/all_samples.vcf (merged)
set -u

REF="/media/desk16/iy16/Downloads/genome/msgenome/sativazm1/zm1chr.fasta"
G=/media/desk16/iy16/Downloads/genome2vcf/target/release/genodiff
OUT=/media/desk16/iy16/Downloads/genome2vcf/batch
THREADS=12

# label=path  (label becomes the sample name in the merged VCF)
QUERIES=(
  "caerulea_GWHJWKE=/media/desk16/iy16/Downloads/genome/msgenome/caerulea/GWHJWKE00000000.1.genome.Chr1-8.fa"
  "caerulea_Mca_long=/media/desk16/iy16/Downloads/genome/msgenome/caerulea/Mca_long.fa"
  "caerulea_mscearula=/media/desk16/iy16/Downloads/genome/msgenome/caerulea/mscearula.fasta"
  "falcata_hap1=/media/desk16/iy16/Downloads/genome/msgenome/falcata/msfalcata_chromosomes.hap1.fa"
  "falcata_hap2=/media/desk16/iy16/Downloads/genome/msgenome/falcata/msfalcata_chromosomes.hap2.fa"
  "falcata_hap3=/media/desk16/iy16/Downloads/genome/msgenome/falcata/msfalcata_chromosomes.hap3.fa"
  "falcata_hap4=/media/desk16/iy16/Downloads/genome/msgenome/falcata/msfalcata_chromosomes.hap4.fa"
  "falcata_GWHHNGT=/media/desk16/iy16/Downloads/genome/msgenome/falcata/GWHHNGT00000000.1.genome.Chr1-8.fa"
  "sativazm1_ZM4_hap1=/media/desk16/iy16/Downloads/genome/msgenome/sativazm1/ZM4.hap1.fa"
  "sativazm1_ZM4_hap2=/media/desk16/iy16/Downloads/genome/msgenome/sativazm1/ZM4.hap2.fa"
  "sativazm1_ZM4_hap3=/media/desk16/iy16/Downloads/genome/msgenome/sativazm1/ZM4.hap3.fa"
  "sativazm1_ZM4_hap4=/media/desk16/iy16/Downloads/genome/msgenome/sativazm1/ZM4.hap4.fa"
  "XJDY_hap1=/media/desk16/iy16/Downloads/genome/XJDY/fasta/XJDY.hap1.fa"
  "XJDY_hap2=/media/desk16/iy16/Downloads/genome/XJDY/fasta/XJDY.hap2.fa"
  "XJDY_hap3=/media/desk16/iy16/Downloads/genome/XJDY/fasta/XJDY.hap3.fa"
  "XJDY_hap4=/media/desk16/iy16/Downloads/genome/XJDY/fasta/XJDY.hap4.fa"
  "zm3_genome=/media/desk16/iy16/Downloads/genome/zm3ver2/fasta/genome.fa"
  "sativazm1_B32=/media/desk16/iy16/Downloads/genome/msgenome/sativazm1/B32.chrom.Chr1-8.fa"
)

mkdir -p "$OUT/links" "$OUT/logs"

# free the previous merged output before starting (it is the largest artifact)
rm -f "$OUT"/all_samples*.vcf

n_ok=0; n_fail=0
t0=$(date +%s)
for entry in "${QUERIES[@]}"; do
  name="${entry%%=*}"
  path="${entry#*=}"
  if [ ! -f "$path" ]; then
    echo "[SKIP] $name: file not found: $path"
    continue
  fi
  ln -sf "$path" "$OUT/links/${name}.fa"
  echo "=== [$name] calling ... $(date +%H:%M:%S)"
  if "$G" call "$REF" "$OUT/links/${name}.fa" -o "$OUT/${name}.vcf" \
       --sample "$name" --threads "$THREADS" > "$OUT/logs/${name}.log" 2>&1; then
    echo "    [ok] $name  vcf=$(du -h "$OUT/${name}.vcf" | cut -f1)  disk: $(df -h "$OUT" | tail -1 | awk '{print $4}') free"
    n_ok=$((n_ok + 1))
  else
    echo "    [FAIL] $name (see $OUT/logs/${name}.log)"
    n_fail=$((n_fail + 1))
  fi
done
echo "batch done in $(( $(date +%s) - t0 ))s: ok=$n_ok fail=$n_fail"

if [ "$n_ok" -lt 2 ]; then
  echo "fewer than 2 VCFs produced, skipping merge"
  exit 1
fi

python3 "$(dirname "$0")/merge_vcfs.py" \
  -o "$OUT/all_samples.vcf" "$OUT"/*.vcf
