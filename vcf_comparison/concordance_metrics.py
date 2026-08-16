#!/usr/bin/env python3
"""
concordance_metrics.py
======================
MIP vs resequencing SNP concordance, reporting the metrics used in
GATK/Picard GenotypeConcordance:

  * genotype confusion matrix (truth = resequencing, call = MIP)
  * sensitivity and PPV for HET and HOM_VAR
  * overall genotype concordance
  * non-reference genotype concordance

Definitions (Picard GenotypeConcordance convention), per genotype class C:
    sensitivity(C) = N(truth=C and call=C) / N(truth=C)
    PPV(C)         = N(truth=C and call=C) / N(call=C)
Only genotypes called on BOTH platforms are counted (missing on either
side is excluded and reported separately).

Genotypes are compared by nucleotide, so a REF/ALT swap between the two
VCFs does not create a false mismatch.

Usage:
    python concordance_metrics.py --mip mip.common.vcf.gz \
        --reseq reseq.common.vcf.gz --out mip_vs_YR_HP #prefix of outputs
"""
import argparse
from collections import defaultdict

import pysam

CLASSES = ["HOM_REF", "HET", "HOM_VAR"]


def is_biallelic_snp(rec):
    if rec.alts is None or len(rec.alts) != 1:
        return False
    return (len(rec.ref) == 1 and len(rec.alts[0]) == 1
            and rec.ref in "ACGT" and rec.alts[0] in "ACGT")


def genotype_nt(rec, sample, min_dp=None):
    call = rec.samples[sample]
    gt = call.get("GT")
    if gt is None or any(a is None for a in gt):
        return None
    if min_dp is not None:
        dp = call.get("DP")
        if dp is None or dp < min_dp:
            return None
    try:
        return tuple(sorted(rec.alleles[a] for a in gt))
    except IndexError:
        return None


def gt_class(gt_nt, ref):
    """Classify a nucleotide genotype relative to the reference allele."""
    n_alt = sum(1 for a in gt_nt if a != ref)
    return CLASSES[n_alt] if n_alt <= 2 else "HOM_VAR"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mip", required=True)
    ap.add_argument("--reseq", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-dp", type=int, default=None)
    args = ap.parse_args()

    mip_vcf = pysam.VariantFile(args.mip)
    reseq_vcf = pysam.VariantFile(args.reseq)
    pairs = [(s, s) for s in mip_vcf.header.samples
             if s in set(reseq_vcf.header.samples)]
    if not pairs:
        raise SystemExit("no shared samples")

    mip_sites = {}
    for rec in mip_vcf:
        if is_biallelic_snp(rec):
            mip_sites[(rec.chrom, rec.pos)] = {
                "pair": frozenset((rec.ref, rec.alts[0])),
                "ref": rec.ref,
                "gts": {m: genotype_nt(rec, m, args.min_dp) for m, _ in pairs},
            }

    matrix = defaultdict(int)          # (truth_class, call_class) -> n
    n_sites = n_allele_mismatch = 0
    missing_mip = missing_reseq = 0
    seen = set()

    for rec in reseq_vcf:
        key = (rec.chrom, rec.pos)
        if key not in mip_sites or not is_biallelic_snp(rec) or key in seen:
            continue
        if frozenset((rec.ref, rec.alts[0])) != mip_sites[key]["pair"]:
            n_allele_mismatch += 1
            continue
        seen.add(key)
        n_sites += 1
        ref = rec.ref
        for m, r in pairs:
            g_mip = mip_sites[key]["gts"][m]
            g_res = genotype_nt(rec, r, args.min_dp)
            if g_res is None:
                missing_reseq += 1
                continue
            if g_mip is None:
                missing_mip += 1
                continue
            matrix[(gt_class(g_res, ref), gt_class(g_mip, ref))] += 1

    total = sum(matrix.values())
    diag = sum(matrix[(c, c)] for c in CLASSES)
    nonref = sum(n for (t, c), n in matrix.items()
                 if not (t == "HOM_REF" and c == "HOM_REF"))
    nonref_ok = sum(matrix[(c, c)] for c in ("HET", "HOM_VAR"))

    def pct(a, b):
        return f"{100.0*a/b:.2f}%" if b else "NA"

    lines = []
    lines.append("Genotype confusion matrix (rows = resequencing truth, cols = MIP call)")
    lines.append("truth\\call\t" + "\t".join(CLASSES) + "\ttruth_total")
    for t in CLASSES:
        row = [matrix[(t, c)] for c in CLASSES]
        lines.append(f"{t}\t" + "\t".join(map(str, row)) + f"\t{sum(row)}")
    lines.append("call_total\t" +
                 "\t".join(str(sum(matrix[(t, c)] for t in CLASSES)) for c in CLASSES) +
                 f"\t{total}")
    lines.append("")
    lines.append("class\tsensitivity\tPPV\tconcordance_n")
    for c in CLASSES:
        n_truth = sum(matrix[(c, x)] for x in CLASSES)
        n_call = sum(matrix[(x, c)] for x in CLASSES)
        both = matrix[(c, c)]
        lines.append(f"{c}\t{pct(both, n_truth)}\t{pct(both, n_call)}\t{both}")
    lines.append("")
    lines.append(f"sites_compared\t{n_sites}")
    lines.append(f"allele_mismatch_sites\t{n_allele_mismatch}")
    lines.append(f"genotypes_compared\t{total}")
    lines.append(f"missing_in_MIP\t{missing_mip}")
    lines.append(f"missing_in_reseq\t{missing_reseq}")
    lines.append(f"overall_genotype_concordance\t{pct(diag, total)}\t({diag}/{total})")
    lines.append(f"non_reference_concordance\t{pct(nonref_ok, nonref)}\t({nonref_ok}/{nonref})")

    text = "\n".join(lines)
    print(text)
    with open(f"{args.out}.metrics.tsv", "w") as fh:
        fh.write(text + "\n")
    print(f"\nwritten: {args.out}.metrics.tsv")


if __name__ == "__main__":
    main()
