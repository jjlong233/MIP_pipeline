#!/usr/bin/env python3
"""
mip_genotype_pipeline.py
========================
End-to-end MIP genotyping pipeline:

    raw paired FASTQ  +  primer file  +  reference genome   ->   formal VCF

This pipeline contains 6 steps in total:

    Stage 0  merge paired-end reads          pandaseq / fastp            [external]
    Stage 1  demultiplex reads by locus      <- 1_seperate_loci_by_primer.py
    Stage 2  group reads by UMI              <- 2_sorted_pe_umi.py
    Stage 3  (optional) top-3 QC report      <- 3_generate_top_three_consensus_seq.py
    Stage 4  per-UMI consensus               <- 4_generate_consensus_sequence.py
    Stage 5  map consensus to reference      bwa mem -> sorted BAM        [external]
    Stage 6  call variants -> VCF            bcftools mpileup + call      [external]

External tools expected on PATH: pandaseq (or fastp), bwa, samtools, bcftools.
Install with e.g.:  conda install -c bioconda pandaseq fastp bwa samtools bcftools

Primer file format (tab-separated, no header):
    <locus_name>\t<primer_L>\t<primer_R>
Both primers are written as ordered (5'->3'). 

Author: Junlong Jiang
"""

import argparse
import glob
import os
import re
import shlex
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

from Bio import SeqIO


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def run(cmd, shell=False):
    """Run an external command, raising on non-zero exit."""
    printable = cmd if isinstance(cmd, str) else " ".join(cmd)
    log(f"  [run] {printable}")
    subprocess.run(cmd, shell=shell, check=True)


def need_tool(name):
    if shutil.which(name) is None:
        sys.exit(
            f"ERROR: required tool '{name}' not found on PATH.\n"
            f"       Install it, e.g.  conda install -c bioconda {name}"
        )


def mkdir(path):
    os.makedirs(path, exist_ok=True)
    return path


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


# --------------------------------------------------------------------------- #
# Sample discovery (batch mode)
# --------------------------------------------------------------------------- #
_FASTQ_EXT = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def _strip_fastq_ext(name):
    for ext in _FASTQ_EXT:
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def discover_samples(sample_dir, r1_token="_R1", r2_token="_R2"):
    """
    Find paired FASTQ files in a folder and pair them into samples.

    A file is treated as R1 if it contains `r1_token`; its mate is the same
    name with r1_token -> r2_token. The sample name is everything before the
    token, with common Illumina suffixes (_S12, _L001, _001) stripped.

    Returns a sorted list of (sample_name, r1_path, r2_path).
    """
    r1_files = []
    for ext in _FASTQ_EXT:
        r1_files += glob.glob(os.path.join(sample_dir, f"*{r1_token}*{ext}"))
    r1_files = sorted(set(r1_files))

    samples, seen = [], {}
    for r1 in r1_files:
        base = os.path.basename(r1)
        r2 = os.path.join(sample_dir, base.replace(r1_token, r2_token, 1))
        if not os.path.exists(r2):
            log(f"  WARNING: no R2 mate for {base} (expected {os.path.basename(r2)}); skipping")
            continue
        stem = _strip_fastq_ext(base)
        name = stem.split(r1_token, 1)[0]
        # strip trailing Illumina tokens (_S12, _L001, _001) in any order
        while True:
            stripped = re.sub(r"_(S\d+|L\d{3}|\d{3})$", "", name)
            if stripped == name:
                break
            name = stripped
        name = name.rstrip("_.") or stem
        if name in seen:
            log(f"  WARNING: duplicate sample name '{name}' "
                f"({base} vs {os.path.basename(seen[name])}); disambiguating")
            name = stem
        seen[name] = r1
        samples.append((name, r1, r2))
    return sorted(samples)


def read_samplesheet(path):
    """
    Explicit sample list. Tab- or comma-separated, optional header. Columns:
        sample_name   r1_path   r2_path
    Returns a list of (sample_name, r1_path, r2_path).
    """
    samples = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[\t,]", line)
            if len(parts) < 3:
                continue
            if parts[0].lower() in ("sample", "sample_name", "name"):
                continue  # header
            samples.append((parts[0].strip(), parts[1].strip(), parts[2].strip()))
    return samples


# --------------------------------------------------------------------------- #
# Stage 0 : merge paired-end reads
# --------------------------------------------------------------------------- #
def merge_reads(r1, r2, out_fasta, merger="pandaseq", threads=4,
                min_len=0, extra_args=""):
    """Merge overlapping read pairs into a single FASTA of merged amplicons."""
    extra = shlex.split(extra_args) if extra_args else []

    if merger == "pandaseq":
        need_tool("pandaseq")
        logfile = out_fasta + ".pandaseq.log"
        # pandaseq writes assembled FASTA to stdout; redirect to file.
        cmd = ["pandaseq", "-f", r1, "-r", r2, "-T", str(threads)]
        if min_len:
            cmd += ["-l", str(min_len)]
        cmd += extra
        full = f"{' '.join(shlex.quote(c) for c in cmd)} > {shlex.quote(out_fasta)} 2> {shlex.quote(logfile)}"
        run(full, shell=True)

    elif merger == "fastp":
        need_tool("fastp")
        merged_fq = out_fasta + ".merged.fastq"
        cmd = ["fastp", "-i", r1, "-I", r2,
               "--merge", "--merged_out", merged_fq,
               "--disable_adapter_trimming",          # keep MIP arms intact
               "-w", str(min(threads, 16)),
               "-j", out_fasta + ".fastp.json",
               "-h", out_fasta + ".fastp.html"]
        if min_len:
            cmd += ["--length_required", str(min_len)]
        cmd += extra
        run(cmd)
        SeqIO.convert(merged_fq, "fastq", out_fasta, "fasta")

    else:
        sys.exit(f"Unknown merger: {merger}")

    n = sum(1 for _ in SeqIO.parse(out_fasta, "fasta"))
    log(f"  merged reads: {n}")
    return n


# --------------------------------------------------------------------------- #
# Stage 1 : demultiplex reads by locus  (from 1_seperate_loci_by_primer.py)
# --------------------------------------------------------------------------- #
def read_primers(primer_file):
    """
    Parse '<locus>\\t<primer_L>\\t<primer_R>' TSV. Both primers are given as
    ordered (5'->3'). primer_L is searched verbatim on the forward/merged strand;
    primer_R is reverse-complemented here so it matches the 3' end of the merged
    read (the reverse oligo's footprint on the forward strand).
    """
    primers = {}
    with open(primer_file) as fh:
        for row in fh:
            if not row.strip():
                continue
            cols = row.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            locus = cols[0]
            pL = cols[1].strip().upper()
            pR = reverse_complement(cols[2].strip().upper())
            primers[locus] = (pL, pR)
    return primers


def iter_fasta(path):
    """Fast minimal FASTA reader (id, sequence). Faster than SeqIO for big files."""
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split(None, 1)[0].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def build_primer_index(primer_list, max_mismatch):
    """
    Build lookup structures for ONE end (all forward primers, or all reverse
    primers). primer_list is [(locus, primer_seq)].

    Returns (exact, lengths, kmer_index, k, by_locus):
      exact       : {primer_seq: [locus, ...]}     -- exact-match fast path
      lengths     : sorted distinct primer lengths
      kmer_index  : {kmer: {locus, ...}}           -- pigeonhole prefilter
      k           : seed length (<= floor(min_len/(max_mismatch+1)) so that at
                    least one seed survives <=max_mismatch substitutions)
      by_locus    : {locus: primer_seq}            -- for fuzzy verification
    """
    exact = defaultdict(list)
    lengths = set()
    by_locus = {}
    min_len = min(len(p) for _, p in primer_list)
    k = max(4, min_len // (max_mismatch + 1))
    kmer_index = defaultdict(set)
    for locus, primer in primer_list:
        exact[primer].append(locus)
        lengths.add(len(primer))
        by_locus[locus] = primer
        if len(primer) < k:
            kmer_index[primer].add(locus)
        else:
            for i in range(len(primer) - k + 1):
                kmer_index[primer[i:i + k]].add(locus)
    return dict(exact), sorted(lengths), dict(kmer_index), k, by_locus


def _hamming_ok(region, start, primer, max_mismatch):
    mism = 0
    for j in range(len(primer)):
        if region[start + j] != primer[j]:
            mism += 1
            if mism > max_mismatch:
                return False
    return True


def match_end(region, index, max_mismatch):
    """Return the set of loci whose primer occurs in `region` (a read end)."""
    exact, lengths, kmer_index, k, by_locus = index
    hits = set()
    R = len(region)

    # 1) exact fast path: slide each distinct primer length, O(1) dict lookups,
    #    cost independent of the number of primers.
    for m in lengths:
        for i in range(R - m + 1):
            loci = exact.get(region[i:i + m])
            if loci:
                hits.update(loci)

    # 2) fuzzy fallback, but only for candidate loci that share a seed k-mer
    #    with this region (pigeonhole: a true <=max_mismatch hit must share one).
    candidates = set()
    for i in range(R - k + 1):
        loci = kmer_index.get(region[i:i + k])
        if loci:
            candidates.update(loci)
    for locus in candidates - hits:
        primer = by_locus[locus]
        m = len(primer)
        for i in range(R - m + 1):
            if _hamming_ok(region, i, primer, max_mismatch):
                hits.add(locus)
                break
    return hits


def demultiplex_by_locus(merged_fasta, primers, max_mismatch=2,
                         skip_N=True, assign="all", umi_len=6, margin=20):
    """
    Assign each merged read to a locus when BOTH primers match. Forward primers
    are searched only in the first `fwd_win` bp and reverse (revcomp) primers in
    the last `rev_win` bp, where the windows cover UMI + primer + margin. This
    keeps behaviour equivalent to the old whole-read scan for real amplicons
    (primers sit at the ends) while being far faster.

    assign='all'   -> a read may land in every matching locus (original behaviour)
    assign='first' -> a read is placed only in the first matching locus
    """
    fwd_list = [(locus, pL) for locus, (pL, _pR) in primers.items()]
    rev_list = [(locus, pR) for locus, (_pL, pR) in primers.items()]
    fwd_index = build_primer_index(fwd_list, max_mismatch)
    rev_index = build_primer_index(rev_list, max_mismatch)

    fwd_win = umi_len + max(fwd_index[1]) + margin
    rev_win = umi_len + max(rev_index[1]) + margin
    order = list(primers.keys())  # preserves insertion order for assign='first'

    buckets = defaultdict(list)
    total = kept = 0
    for rid, seq in iter_fasta(merged_fasta):
        total += 1
        seq = seq.upper()
        if skip_N and "N" in seq:
            continue

        fwd_hits = match_end(seq[:fwd_win], fwd_index, max_mismatch)
        if not fwd_hits:                      # no forward primer -> unassignable
            continue
        rev_hits = match_end(seq[-rev_win:], rev_index, max_mismatch)
        both = fwd_hits & rev_hits
        if not both:
            continue

        if assign == "first":
            for locus in order:
                if locus in both:
                    buckets[locus].append((rid, seq))
                    break
        else:
            for locus in both:
                buckets[locus].append((rid, seq))
        kept += 1

    log(f"  reads assigned to a locus: {kept}/{total}")
    for locus in sorted(buckets):
        log(f"    {locus}: {len(buckets[locus])} reads")
    return buckets


# --------------------------------------------------------------------------- #
# Stages 2 + 4 fused : UMI grouping + per-UMI consensus
#   (from 2_sorted_pe_umi.py and 4_generate_consensus_sequence.py)
# --------------------------------------------------------------------------- #
def consensus_by_umi(reads, umi_len=6, min_reads=3, trim=6):
    """
    reads : [(read_id, sequence)] for ONE locus.

    Group by (first umi_len bp, last umi_len bp). For each UMI family whose
    total read count >= min_reads, take the most abundant *exact* sequence as
    that molecule's consensus, then trim `trim` bp from each end (removes UMIs).

    Returns
        consensus : [(umi_f, umi_r, family_size, consensus_seq)]
        families  : {(umi_f, umi_r): Counter(seq -> count)}   (for the QC report)
    """
    families = defaultdict(Counter)
    for _rid, seq in reads:
        umi = (seq[:umi_len], seq[-umi_len:])
        families[umi][seq] += 1

    consensus = []
    for (uf, ur), counter in families.items():
        family_size = sum(counter.values())
        if family_size < min_reads:
            continue
        top_seq = counter.most_common(1)[0][0]
        trimmed = top_seq[trim:-trim] if trim > 0 else top_seq
        if trimmed:
            consensus.append((uf, ur, family_size, trimmed))
    return consensus, families


# --------------------------------------------------------------------------- #
# Stage 3 : optional QC report  (from 3_generate_top_three_consensus_seq.py)
# --------------------------------------------------------------------------- #
def write_top_three_report(families, out_path, top=3):
    total = sum(sum(c.values()) for c in families.values()) or 1
    with open(out_path, "w") as fh:
        for (uf, ur), counter in families.items():
            fam = sum(counter.values())
            fh.write(f"UMI_F: {uf}\n")
            fh.write(f"UMI_R: {ur}\n")
            fh.write(f"Frequency: {fam} (Frequency Percentage: {fam / total:.2%})\n")
            fh.write("Sequences:\n")
            for seq, cnt in counter.most_common(top):
                fh.write(f"{seq} (Read Count: {cnt}) (Percentage: {cnt / fam:.2%})\n")
            fh.write("\n")


# --------------------------------------------------------------------------- #
# Consensus writers
# --------------------------------------------------------------------------- #
def write_consensus(all_consensus, sample, fasta_path, fastq_path, qual="I"):
    """
    all_consensus : [(locus, umi_f, umi_r, family_size, seq)]

    Writes BOTH:
      * a FASTA (your workflow's deliverable), and
      * a FASTQ used for mapping. bcftools needs per-base qualities, so we give
        every consensus base a fixed high quality ('I' == Q40). A consensus base
        drawn from a UMI family is high-confidence, so a uniform high Q is a
        reasonable stand-in and keeps the pileup step clean.
    Read names are made unique + space-free so BWA/SAM accept them.
    """
    with open(fasta_path, "w") as fa, open(fastq_path, "w") as fq:
        for locus, uf, ur, size, seq in all_consensus:
            name = f"{sample}|{locus}|{uf}_{ur}|n{size}"
            fa.write(f">{name}\n{seq}\n")
            fq.write(f"@{name}\n{seq}\n+\n{qual * len(seq)}\n")


# --------------------------------------------------------------------------- #
# Stage 5 : map consensus to reference
# --------------------------------------------------------------------------- #
def ensure_reference_index(reference):
    need_tool("bwa")
    need_tool("samtools")
    if not os.path.exists(reference + ".bwt"):
        log("  indexing reference with bwa (one-time)...")
        run(["bwa", "index", reference])
    if not os.path.exists(reference + ".fai"):
        run(["samtools", "faidx", reference])


def map_consensus(consensus_fastq, reference, out_bam, sample, threads=4):
    rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:{sample}_MIP"
    cmd = (
        f"bwa mem -t {threads} -R '{rg}' "
        f"{shlex.quote(reference)} {shlex.quote(consensus_fastq)} "
        f"| samtools sort -@ {threads} -o {shlex.quote(out_bam)} -"
    )
    run(cmd, shell=True)
    run(["samtools", "index", out_bam])


# --------------------------------------------------------------------------- #
# Stage 6 : variant calling -> VCF
# --------------------------------------------------------------------------- #
def call_variants(bams, reference, out_prefix, threads=4,
                  max_depth=1000000, split=True):
    """
    One or more sorted BAMs -> bcftools mpileup | bcftools call -> normalized VCF.

    Passing several BAMs (each carrying its own @RG SM tag) makes bcftools emit
    a single JOINT multi-sample VCF with one genotype column per sample. This
    replaces hand-parsing CIGAR/MD:Z: bcftools produces a spec-compliant VCF with
    QUAL, per-sample GT, DP and per-allele AD, left-aligns indels, and infers
    homozygote/heterozygote from the allele balance across UMI-consensus
    molecules -- folding your 'consensus_filter' het/hom step into the caller.
    """
    need_tool("bcftools")
    if isinstance(bams, str):
        bams = [bams]
    raw = f"{out_prefix}.raw.vcf.gz"
    norm = f"{out_prefix}.norm.vcf.gz"
    final_vcf = f"{out_prefix}.vcf"

    # a bam list keeps the command short and identical for 1 or N samples
    bamlist = f"{out_prefix}.bam.list"
    with open(bamlist, "w") as fh:
        for b in bams:
            fh.write(b + "\n")

    mpileup = (
        f"bcftools mpileup -f {shlex.quote(reference)} -a AD,DP "
        f"--max-depth {max_depth} --threads {threads} -b {shlex.quote(bamlist)}"
    )
    call = f"bcftools call -m -v --threads {threads} -Oz -o {shlex.quote(raw)}"
    run(f"{mpileup} | {call}", shell=True)
    run(["bcftools", "index", raw])

    # left-align indels + split multiallelic records into biallelic
    run(["bcftools", "norm", "-f", reference, "-m", "-both",
         "-Oz", "-o", norm, raw])
    run(["bcftools", "index", norm])

    # plain-text VCF deliverable
    run(["bcftools", "view", norm, "-o", final_vcf])

    outputs = {"vcf": final_vcf, "vcf_gz": norm}
    if split:
        snps = f"{out_prefix}.snps.vcf"
        indels = f"{out_prefix}.indels.vcf"
        run(["bcftools", "view", "-v", "snps", norm, "-o", snps])
        run(["bcftools", "view", "-v", "indels", norm, "-o", indels])
        outputs["snps"] = snps
        outputs["indels"] = indels

    n_var = 0
    with open(final_vcf) as fh:
        for line in fh:
            if not line.startswith("#"):
                n_var += 1
    log(f"  variants called: {n_var}  (samples: {len(bams)})")
    return outputs


# --------------------------------------------------------------------------- #
# Per-sample processing: stages 0-5  (returns the sorted BAM, or None)
# --------------------------------------------------------------------------- #
def process_sample(sample, r1, r2, merged, primers, reference, sample_outdir, args):
    """Run merge -> demux -> UMI consensus -> map for ONE sample.

    Returns (bam_path, consensus_fasta) on success, or None if the sample yielded
    no consensus molecules (so it can be skipped in a batch without aborting)."""
    d_merge = mkdir(os.path.join(sample_outdir, "01.merged"))
    d_loci = mkdir(os.path.join(sample_outdir, "02.loci"))
    d_qc = mkdir(os.path.join(sample_outdir, "03.qc"))
    d_cons = mkdir(os.path.join(sample_outdir, "04.consensus"))
    d_map = mkdir(os.path.join(sample_outdir, "05.mapping"))

    log(f"--- sample: {sample} ---")

    # Stage 0 : merge
    if merged:
        merged_fasta = merged
        log("  Stage 0: using provided merged FASTA")
    else:
        merged_fasta = os.path.join(d_merge, f"{sample}.merged.fasta")
        log(f"  Stage 0: merging reads with {args.merger}")
        merge_reads(r1, r2, merged_fasta, merger=args.merger,
                    threads=args.threads, min_len=args.min_len,
                    extra_args=args.merger_args)

    # Stage 1 : demultiplex
    log("  Stage 1: demultiplexing by locus")
    buckets = demultiplex_by_locus(
        merged_fasta, primers,
        max_mismatch=args.max_mismatch,
        skip_N=not args.keep_N,
        assign=args.assign,
        umi_len=args.umi_len,
    )

    # Stages 2-4 : UMI grouping + consensus
    log("  Stages 2-4: UMI grouping + consensus")
    all_consensus = []
    summary_rows = []
    for locus, reads in buckets.items():
        with open(os.path.join(d_loci, f"{locus}.fasta"), "w") as fh:
            for rid, seq in reads:
                fh.write(f">{rid}\n{seq}\n")
        consensus, families = consensus_by_umi(
            reads, umi_len=args.umi_len,
            min_reads=args.min_reads, trim=args.trim,
        )
        for uf, ur, size, seq in consensus:
            all_consensus.append((locus, uf, ur, size, seq))
        summary_rows.append((locus, len(reads), len(families), len(consensus)))
        if args.qc_report:
            write_top_three_report(
                families, os.path.join(d_qc, f"{locus}_top_three.txt"))

    with open(os.path.join(d_cons, f"{sample}.consensus_summary.tsv"), "w") as fh:
        fh.write("locus\treads\tumi_families\tconsensus_molecules\n")
        for row in sorted(summary_rows):
            fh.write("\t".join(map(str, row)) + "\n")
    log(f"  total consensus molecules: {len(all_consensus)}")

    if not all_consensus:
        log(f"  WARNING: sample '{sample}' produced no consensus; skipping it.")
        return None

    fasta_path = os.path.join(d_cons, f"{sample}.consensus.fasta")
    fastq_path = os.path.join(d_cons, f"{sample}.consensus.fastq")
    write_consensus(all_consensus, sample, fasta_path, fastq_path)

    # Stage 5 : map
    log("  Stage 5: mapping consensus to reference")
    bam = os.path.join(d_map, f"{sample}.consensus.sorted.bam")
    map_consensus(fastq_path, reference, bam, sample, threads=args.threads)
    return bam, fasta_path


def main():
    ap = argparse.ArgumentParser(
        description="MIP genotyping: raw paired FASTQ + primers + reference -> VCF "
                    "(single sample, or batch a folder into one joint VCF)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # --- input: choose ONE of these three modes ---
    ap.add_argument("--r1", help="[single] Raw forward FASTQ")
    ap.add_argument("--r2", help="[single] Raw reverse FASTQ")
    ap.add_argument("--merged", help="[single] Skip merging: already-merged FASTA")
    ap.add_argument("--sample-dir", help="[batch] Folder of paired FASTQs; auto-pairs samples")
    ap.add_argument("--samplesheet", help="[batch] TSV/CSV: sample_name<TAB>r1<TAB>r2")
    ap.add_argument("--r1-token", default="_R1", help="[batch] token marking R1 files")
    ap.add_argument("--r2-token", default="_R2", help="[batch] token marking R2 files")
    # --- shared required inputs ---
    ap.add_argument("--primer", required=True,
                    help="Primer TSV: locus<TAB>primer_L<TAB>primer_R, both 5'->3' as ordered")
    ap.add_argument("--reference", required=True, help="Reference genome FASTA")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--sample", help="[single] Sample name (default: from --r1/--merged)")
    # stage 0
    ap.add_argument("--merger", choices=["pandaseq", "fastp"], default="pandaseq")
    ap.add_argument("--merger-args", default="", help="Extra args passed to the merger")
    ap.add_argument("--min-len", type=int, default=0, help="Minimum merged length")
    # stages 1-4
    ap.add_argument("--umi-len", type=int, default=6, help="UMI length at each end")
    ap.add_argument("--min-reads", type=int, default=3, help="Min reads per UMI family")
    ap.add_argument("--trim", type=int, default=6, help="bp trimmed from each end (UMI removal)")
    ap.add_argument("--max-mismatch", type=int, default=2, help="Allowed primer mismatches")
    ap.add_argument("--assign", choices=["all", "first"], default="all",
                    help="Assign a read to all matching loci or only the first")
    ap.add_argument("--keep-N", action="store_true", help="Keep reads containing N")
    ap.add_argument("--qc-report", action="store_true", help="Write stage-3 top-3 QC reports")
    # stages 5-6
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--no-split", action="store_true", help="Do not split SNP/indel VCFs")
    args = ap.parse_args()

    batch = bool(args.sample_dir or args.samplesheet)
    if not batch and not args.merged and not (args.r1 and args.r2):
        ap.error("Provide a single sample (--r1/--r2 or --merged) "
                 "OR a batch (--sample-dir or --samplesheet).")

    out = mkdir(args.outdir)
    primers = read_primers(args.primer)
    log(f"loaded {len(primers)} loci from {args.primer}")
    ensure_reference_index(args.reference)   # once, shared by all samples

    # ---- build the list of samples to process ----
    if args.samplesheet:
        sample_list = read_samplesheet(args.samplesheet)
    elif args.sample_dir:
        sample_list = discover_samples(args.sample_dir, args.r1_token, args.r2_token)
    else:
        sample = args.sample or os.path.basename(args.r1 or args.merged).split(".")[0]
        sample_list = [(sample, args.r1, args.r2)]

    if not sample_list:
        sys.exit("No samples found. Check --sample-dir / --samplesheet / tokens.")

    log(f"=== MIP pipeline | {len(sample_list)} sample(s) | "
        f"{'BATCH -> joint VCF' if batch else 'single'} ===")
    for name, r1, _r2 in sample_list:
        log(f"  - {name}  ({os.path.basename(r1) if r1 else args.merged})")

    # ---- process each sample one by one (stages 0-5) ----
    bams, processed, failed = [], [], []
    for name, r1, r2 in sample_list:
        sample_outdir = os.path.join(out, "samples", name) if batch else out
        mkdir(sample_outdir)
        merged = args.merged if (not batch) else None
        try:
            result = process_sample(name, r1, r2, merged, primers,
                                    args.reference, sample_outdir, args)
        except subprocess.CalledProcessError as e:
            log(f"  ERROR processing '{name}': {e}. Skipping.")
            failed.append(name)
            continue
        if result is None:
            failed.append(name)
            continue
        bam, _fasta = result
        bams.append(bam)
        processed.append(name)

    if not bams:
        sys.exit("No sample produced a BAM; nothing to call. Check inputs / --min-reads.")

    # ---- joint variant calling (stage 6) ----
    log(f"Stage 6: {'joint ' if len(bams) > 1 else ''}variant calling "
        f"over {len(bams)} sample(s)")
    vcf_dir = mkdir(os.path.join(out, "joint")) if batch else mkdir(os.path.join(out, "06.vcf"))
    prefix = os.path.join(vcf_dir, "cohort" if batch else processed[0])
    outputs = call_variants(
        bams, args.reference, prefix,
        threads=args.threads, split=not args.no_split,
    )

    log("=== done ===")
    log(f"Samples processed : {len(processed)} ({', '.join(processed)})")
    if failed:
        log(f"Samples skipped   : {len(failed)} ({', '.join(failed)})")
    log(f"Joint VCF         : {outputs['vcf']}"
        if batch else f"Final VCF         : {outputs['vcf']}")
    if "snps" in outputs:
        log(f"SNP VCF           : {outputs['snps']}")
        log(f"Indel VCF         : {outputs['indels']}")


if __name__ == "__main__":
    main()
