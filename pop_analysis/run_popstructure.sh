#!/usr/bin/env bash
# =============================================================================
# run_popstructure.sh
# This script is for generating "popstructure" folder needed by "population_analysis_plot.py"
#
# Follows the manuscript Methods exactly:
#   VCFtools v0.1.16 : biallelic SNPs, MAF >= 5%, per-site missingness <= 20%
#   PLINK v1.9       : binary files; LD prune 50-SNP window, 10-SNP step, r2 > 0.5
#   PCA              : from the LD-pruned set (eigenvalues -> % variance)
#   PLINK            : remove HWE p < 0.001 for admixture
#   ADMIXTURE v1.3.0 : K = 1-10 with cross-validation
#
# Usage:  bash run_popstructure.sh cohort_renamed.ontarget.vcf.gz pop_map.txt outdir
# =============================================================================
set -euo pipefail

VCF=${1:-cohort_renamed.ontarget.vcf.gz}
POPMAP=${2:-pop_map.txt}          # 2 cols: sample <TAB> population
OUT=${3:-popstructure}
THREADS=${THREADS:-4}

mkdir -p "$OUT"; cd "$OUT"
VCF=$(readlink -f "../$VCF" 2>/dev/null || readlink -f "$VCF")
POPMAP=$(readlink -f "../$POPMAP" 2>/dev/null || readlink -f "$POPMAP")

echo "### 0. input summary"
bcftools query -l "$VCF" | wc -l | xargs echo "  samples:"
bcftools view -H "$VCF" | wc -l | xargs echo "  sites  :"

# -----------------------------------------------------------------------------
# 1. VCFtools filtering: biallelic SNPs, MAF >= 0.05, missingness <= 20%
#    NOTE: vcftools --max-missing 0.8 KEEPS sites with >=80% call rate,
#          i.e. <=20% missing (0 = allow all missing, 1 = allow none).
# -----------------------------------------------------------------------------
echo "### 1. VCFtools filtering (MAF>=0.05, missingness<=20%, biallelic SNPs)"
vcftools --gzvcf "$VCF" \
    --remove-indels --min-alleles 2 --max-alleles 2 \
    --maf 0.05 --max-missing 0.8 \
    --recode --recode-INFO-all --out mip.filtered
mv mip.filtered.recode.vcf mip.filtered.vcf
echo -n "  SNPs after filtering: "; grep -vc '^#' mip.filtered.vcf

# -----------------------------------------------------------------------------
# 2. PLINK binary. --allow-extra-chr is REQUIRED: contigs are NC_xxxxxx.1,
#    not human chromosome codes. --double-id keeps VCF sample names as FID+IID.
# -----------------------------------------------------------------------------
echo "### 2. PLINK binary"
plink --vcf mip.filtered.vcf --allow-extra-chr --double-id \
      --set-missing-var-ids '@:#' --make-bed --out mip

# -----------------------------------------------------------------------------
# 3. LD pruning: 50-SNP window, 10-SNP step, r2 > 0.5
# -----------------------------------------------------------------------------
echo "### 3. LD pruning (50 10 0.5)"
plink --bfile mip --allow-extra-chr --indep-pairwise 50 10 0.5 --out mip.ld
plink --bfile mip --allow-extra-chr --extract mip.ld.prune.in --make-bed --out mip.pruned
echo -n "  SNPs after LD pruning: "; wc -l < mip.pruned.bim

# -----------------------------------------------------------------------------
# 4. PCA from the pruned set
# -----------------------------------------------------------------------------
echo "### 4. PCA"
plink --bfile mip.pruned --allow-extra-chr --pca 10 --out mip.pca
echo "  -> mip.pca.eigenvec / mip.pca.eigenval"

# -----------------------------------------------------------------------------
# 5. Admixture input: additionally remove HWE p < 0.001
# -----------------------------------------------------------------------------
echo "### 5. HWE filter for admixture (p < 0.001)"
plink --bfile mip.pruned --allow-extra-chr --hwe 0.001 --make-bed --out mip.admix
echo -n "  SNPs after HWE filter: "; wc -l < mip.admix.bim

# ADMIXTURE cannot parse non-numeric chromosome names (NC_066509.1).
# It ignores chromosome information, so recode column 1 to a dummy integer.
cp mip.admix.bim mip.admix.bim.orig
awk 'BEGIN{OFS="\t"}{$1=1; print}' mip.admix.bim.orig > mip.admix.bim

# -----------------------------------------------------------------------------
# 6. ADMIXTURE K = 1..10 with cross-validation
# -----------------------------------------------------------------------------
echo "### 6. ADMIXTURE K=1..10 (cross-validation)"
for K in $(seq 1 10); do
    admixture --cv=10 -j${THREADS} mip.admix.bed $K | tee "admixture_K${K}.log"
done

grep -h "CV error" admixture_K*.log \
  | sed 's/.*(K=\([0-9]*\)).*: \(.*\)/\1\t\2/' | sort -n > cv_error.tsv
echo "### cross-validation error by K"
awk 'BEGIN{print "K\tCV_error"}{print}' cv_error.tsv
BESTK=$(sort -k2,2g cv_error.tsv | head -1 | cut -f1)
echo "### best K (lowest CV error): $BESTK"

# sample order for the Q matrices (ADMIXTURE Q rows follow the .fam order)
cut -d' ' -f2 mip.admix.fam > sample_order.txt
cp "$POPMAP" pop_map.txt

echo
echo "=== DONE ==="
echo "PCA        : $OUT/mip.pca.eigenvec , mip.pca.eigenval"
echo "ADMIXTURE  : $OUT/mip.admix.<K>.Q  (K=1..10)"
echo "CV errors  : $OUT/cv_error.tsv  (best K = $BESTK)"
echo "Next: python plot_fig6.py --dir $OUT --popmap pop_map.txt --out Fig6_new"
