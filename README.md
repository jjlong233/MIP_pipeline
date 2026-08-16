# MIP genotyping pipeline for Chinese mitten crab (*Eriocheir sinensis*)

Scripts accompanying:

> Jiang, J., et al. Development of a high-throughput molecular inversion probe panel for SNP genotyping in Chinese mitten crab (*Eriocheir sinensis*).

A 48-probe molecular inversion probe (MIP) panel and an end-to-end pipeline that converts raw enrichment reads into a normalized multi-sample VCF, without genome-wide variant calling.

---

## Contents

"MIP_genotyping_pipeline.py": Raw FASTQ → multi-sample VCF. Merging, per-locus demultiplexing, UMI consensus, alignment, joint variant calling.
"concordance_metrics.py": Benchmarks MIP genotypes against a whole-genome resequencing call set: concordance, sensitivity and PPV.
"run_popstructure.sh": Filtering, LD pruning, PCA and ADMIXTURE from an on-target VCF.
"population_analysis_plot.py": Publication figure: PCA, admixture bar plots and cross-validation error.

---

## Requirements

**Python** ≥ 3.9 with `numpy`, `pysam`, `matplotlib`

**External tools** on `$PATH`:

| Tool | Version used |
|---|---|
| PANDAseq | 2.11 |
| BWA | 0.7.17 |
| SAMtools | 1.21 |
| BCFtools | 1.5 |
| VCFtools | 0.1.17 |
| PLINK | 1.9 |
| ADMIXTURE | 1.3.0 |

```bash
conda create -n mip -c bioconda -c conda-forge \
    python=3.11 numpy pysam matplotlib \
    pandaseq bwa samtools bcftools vcftools plink admixture
conda activate mip
```

---

## 1. Genotyping pipeline

### Input

A tab-separated probe file, one locus per row, both arms written 5'→3'
(the reverse arm is reverse-complemented internally):

```
locus_id            primer_F                primer_R
LG01_1592493        CACACACACAAAGCTCTCCC    GTTCGTGTTTTGCCTGACCT
LG02_1008402        TTTCGCCTTATCTGCTCCCT    AGTCTTGTGCAGTGGGGTAA
```

### Single sample

```bash
python MIP_genotyping_pipeline.py \
    --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
    --primers probes.tsv \
    --ref GCF_024679095.1_genomic.fna \
    --sample SAMPLE01 \
    --outdir results/SAMPLE01
```

### Batch, with joint calling across the cohort

```bash
python MIP_genotyping_pipeline.py \
    --sample-dir fastq/ \
    --primers probes.tsv \
    --ref GCF_024679095.1_genomic.fna \
    --outdir results \
    --threads 8
```

Or with an explicit sample sheet (`sample_id`, `r1`, `r2`):

```bash
python MIP_genotyping_pipeline.py --samplesheet samples.tsv ... 
```

### Key parameters

| Option | Default | Meaning |
|---|---|---|
| `--umi-len` | 6 | UMI length at each read end |
| `--min-reads` | 3 | Reads required to call a UMI family |
| `--max-mismatch` | 2 | Mismatches tolerated per probe arm |
| `--merger` | pandaseq | Read merger (`pandaseq` or `fastp`) |

Raising `--min-reads` reduces allelic dropout at the cost of missing genotypes:
a heterozygous locus represented by *n* independent molecules yields molecules of
a single allele with probability 2 × 0.5<sup>n</sup> (25% at n = 3, <1% at n = 8).

### Output

```
results/
├── samples/<sample>/
│   ├── 01.merged/          merged reads
│   ├── 02.demux/           per-locus read groups
│   ├── 03.umi/             UMI families
│   ├── 04.consensus/       consensus sequences + consensus_summary.tsv
│   └── 05.align/           sorted, indexed BAM
└── joint/
    ├── cohort.raw.vcf.gz
    ├── cohort.norm.vcf.gz
    ├── cohort.snp.vcf.gz
    └── cohort.indel.vcf.gz
```

Because each consensus sequence derives from one UMI-tagged molecule,
`FORMAT/DP` is the number of independent molecules observed at a site.

---

## 2. Concordance against resequencing

```bash
python concordance_metrics.py \
    --mip cohort.snp.vcf.gz \
    --reseq resequencing.vcf.gz \
    --out mip_vs_reseq
```

A site is compared only when both platforms call a biallelic SNP at the same
position with the same pair of alleles. Genotypes are compared as allele pairs
rather than allele indices, so differences in reference/alternate designation
between call sets do not generate spurious mismatches. Genotypes missing in
either dataset are excluded and tallied separately.

Writes `<out>.summary.txt` and `<out>.metrics.tsv` containing the genotype
confusion matrix, overall and non-reference concordance, and per-class
sensitivity and positive predictive value.

---

## 3. Population structure

```bash
bash run_popstructure.sh cohort.snp.vcf.gz pop_map.txt popstructure
```

`pop_map.txt` is two columns, `sample_id` and `population`:

```
E1401    YR
H01      HP
```

Applies MAF ≥ 0.05 and per-site missingness ≤ 0.2, removes duplicate positions,
prunes for linkage disequilibrium (50-SNP window, 10-SNP step, r² > 0.5), and
runs PCA plus ADMIXTURE for K = 1–10 with cross-validation.

```bash
python population_analysis_plot.py \
    --dir popstructure --popmap pop_map.txt \
    --out Fig_popstructure --cv-panel
```
