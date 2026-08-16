#!/usr/bin/env python3
"""
concordance_mip_vs_reseq.py
===========================
Compare MIP genotypes against (merged) resequencing genotypes at SNPs that
BOTH platforms call, to estimate MIP genotyping accuracy.

Inputs
    --mip     MIP VCF (e.g. cohort.vcf / cohort.norm.vcf.gz)
    --reseq   resequencing VCF, already MERGED across the two files with
              `bcftools merge` (see the header notes / accompanying commands)
    --out     output prefix
    --sample-map  optional TSV mapping MIP sample name <TAB> reseq sample name
                  (use when the same individual is named differently)
    --min-dp  optional: treat a genotype with FORMAT/DP below this as missing
    --min-gq  optional: treat a genotype with FORMAT/GQ below this as missing

What counts as an "overlapping SNP"
    Same CHROM, same POS, both biallelic SNPs, and the unordered allele pair
    {REF, ALT} is identical (nucleotide-level, so a REF/ALT swap between the two
    files still matches). Sites where the two files disagree on the alt allele
    are reported separately as allele mismatches and are NOT scored.

How genotypes are compared
    Each genotype is decoded to its nucleotides using that record's own alleles,
    e.g. reseq 1/0 with REF=G,ALT=A -> ('A','G'); MIP 0/1 with REF=A,ALT=G ->
    ('A','G'). Zygosity is preserved (0/0 -> ('A','A'), 1/1 -> ('G','G')), so
    hom/het mismatches are caught, but REF/ALT order is irrelevant.

Metrics reported
    * overlapping SNP sites, and how many samples were comparable
    * overall genotype concordance = matching / compared genotypes
    * non-reference concordance = same, excluding calls that are hom-ref in both
      (hom-ref/hom-ref dominates and otherwise inflates the number)
    * per-sample and per-site tables, plus every discordant call

Requires: pysam  (pip install pysam)
"""

import argparse
import sys
from collections import defaultdict

import pysam


def load_sample_map(path):
    m = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            a, b = line.split()[:2]
            m[a] = b
    return m


def is_biallelic_snp(rec):
    if rec.alts is None or len(rec.alts) != 1:
        return False
    ref, alt = rec.ref, rec.alts[0]
    return (len(ref) == 1 and len(alt) == 1
            and ref in "ACGT" and alt in "ACGT")


def genotype_nt(rec, sample, min_dp=None, min_gq=None):
    """Return sorted tuple of allele nucleotides, or None if missing/filtered."""
    call = rec.samples[sample]
    gt = call.get("GT")
    if gt is None or any(a is None for a in gt):
        return None
    if min_dp is not None:
        dp = call.get("DP")
        if dp is None or dp < min_dp:
            return None
    if min_gq is not None:
        gq = call.get("GQ")
        if gq is None or gq < min_gq:
            return None
    try:
        return tuple(sorted(rec.alleles[a] for a in gt))
    except IndexError:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mip", required=True)
    ap.add_argument("--reseq", required=True, help="already bcftools-merged reseq VCF")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--sample-map", help="TSV: mip_sample <TAB> reseq_sample")
    ap.add_argument("--min-dp", type=int, default=None)
    ap.add_argument("--min-gq", type=int, default=None)
    args = ap.parse_args()

    mip_vcf = pysam.VariantFile(args.mip)
    reseq_vcf = pysam.VariantFile(args.reseq)
    mip_samples = list(mip_vcf.header.samples)
    reseq_samples = set(reseq_vcf.header.samples)

    # ---- resolve the MIP<->reseq sample correspondence ----
    if args.sample_map:
        smap = load_sample_map(args.sample_map)
        pairs = [(m, r) for m, r in smap.items()
                 if m in mip_samples and r in reseq_samples]
    else:
        pairs = [(s, s) for s in mip_samples if s in reseq_samples]

    if not pairs:
        sys.exit("ERROR: no shared samples between MIP and reseq.\n"
                 "       Supply --sample-map, or run `bcftools gtcheck` to find\n"
                 "       which reseq sample each MIP sample corresponds to.")
    sys.stderr.write(f"comparing {len(pairs)} sample pair(s): "
                     f"{', '.join(m+'~'+r for m, r in pairs[:10])}"
                     f"{' ...' if len(pairs) > 10 else ''}\n")

    # ---- index MIP biallelic SNPs: (chrom,pos) -> {pair, ref, alt, gts} ----
    mip_sites = {}
    for rec in mip_vcf:
        if not is_biallelic_snp(rec):
            continue
        key = (rec.chrom, rec.pos)
        mip_sites[key] = {
            "pair": frozenset((rec.ref, rec.alts[0])),
            "gts": {m: genotype_nt(rec, m, args.min_dp, args.min_gq)
                    for m, _ in pairs},
        }
    sys.stderr.write(f"MIP biallelic SNP sites: {len(mip_sites)}\n")

    # ---- stream reseq, compare at shared positions ----
    n_overlap = 0
    n_allele_mismatch = 0
    n_pos_in_reseq = 0
    compared = concordant = 0
    nonref_compared = nonref_concordant = 0
    per_sample = defaultdict(lambda: [0, 0])        # sample -> [compared, concordant]
    per_site = {}                                    # key -> [compared, concordant]
    discordant = []                                  # (chrom,pos,sample,mip_gt,reseq_gt)
    seen_pos = set()

    for rec in reseq_vcf:
        key = (rec.chrom, rec.pos)
        if key not in mip_sites or not is_biallelic_snp(rec):
            continue
        n_pos_in_reseq += 1
        reseq_pair = frozenset((rec.ref, rec.alts[0]))
        if reseq_pair != mip_sites[key]["pair"]:
            # same position, different SNP allele -> not the same variant
            if key not in seen_pos:
                n_allele_mismatch += 1
            continue
        if key in seen_pos:
            continue                                 # already scored this position
        seen_pos.add(key)
        n_overlap += 1
        per_site[key] = [0, 0]

        for m, r in pairs:
            g_mip = mip_sites[key]["gts"][m]
            g_reseq = genotype_nt(rec, r, args.min_dp, args.min_gq)
            if g_mip is None or g_reseq is None:
                continue
            compared += 1
            per_sample[m][0] += 1
            per_site[key][0] += 1
            match = (g_mip == g_reseq)
            ref_nt = rec.ref
            both_homref = (set(g_mip) == {ref_nt} and set(g_reseq) == {ref_nt})
            if not both_homref:
                nonref_compared += 1
            if match:
                concordant += 1
                per_sample[m][1] += 1
                per_site[key][1] += 1
                if not both_homref:
                    nonref_concordant += 1
            else:
                discordant.append((rec.chrom, rec.pos, m,
                                   "/".join(g_mip), "/".join(g_reseq)))

    # ---- write outputs ----
    def pct(a, b):
        return f"{100.0 * a / b:.2f}%" if b else "NA"

    with open(f"{args.out}.per_sample.tsv", "w") as fh:
        fh.write("sample\tcompared\tconcordant\tconcordance\n")
        for s in sorted(per_sample):
            c, k = per_sample[s]
            fh.write(f"{s}\t{c}\t{k}\t{pct(k, c)}\n")

    with open(f"{args.out}.per_site.tsv", "w") as fh:
        fh.write("chrom\tpos\tcompared\tconcordant\tconcordance\n")
        for (chrom, pos) in sorted(per_site):
            c, k = per_site[(chrom, pos)]
            fh.write(f"{chrom}\t{pos}\t{c}\t{k}\t{pct(k, c)}\n")

    with open(f"{args.out}.discordant.tsv", "w") as fh:
        fh.write("chrom\tpos\tsample\tmip_gt\treseq_gt\n")
        for row in discordant:
            fh.write("\t".join(map(str, row)) + "\n")

    # ---- summary ----
    print("=" * 60)
    print("MIP vs resequencing SNP concordance")
    print("=" * 60)
    print(f"sample pairs compared          : {len(pairs)}")
    print(f"MIP biallelic SNP sites        : {len(mip_sites)}")
    print(f"overlapping SNP sites (scored) : {n_overlap}")
    print(f"same-position allele mismatches: {n_allele_mismatch} (not scored)")
    print(f"genotypes compared             : {compared}")
    print(f"overall genotype concordance   : {pct(concordant, compared)} "
          f"({concordant}/{compared})")
    print(f"non-reference concordance      : {pct(nonref_concordant, nonref_compared)} "
          f"({nonref_concordant}/{nonref_compared})")
    print(f"discordant genotype calls      : {len(discordant)}")
    print("-" * 60)
    print(f"tables written: {args.out}.per_sample.tsv, "
          f"{args.out}.per_site.tsv, {args.out}.discordant.tsv")


if __name__ == "__main__":
    main()
