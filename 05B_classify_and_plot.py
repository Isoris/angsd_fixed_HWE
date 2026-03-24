#!/usr/bin/env python3
"""
05B_classify_and_plot.py

Combined classifier + QC plotter for the ANGSD SNP failure-mode debugging workflow.

PART 7: Classification & BED outputs
PART 8: QC plots from .mafs.gz data

Produces:
  BED tracks:
    - PASS_all_filters.bed         (sites in Pass C)
    - FAIL_not_callable.bed        (outside callable mask — interval BED)
    - FAIL_upstream.bed            (callable but not in Pass A — interval BED)
    - FAIL_minMaf_only.bed         (in Pass B, not Pass C)
    - FAIL_snpPval_only.bed        (in Pass A, not Pass B, MAF >= threshold)
    - FAIL_multiple_filters.bed    (in Pass A, not Pass B, MAF < threshold)
    - A_not_B.bed                  (Pass A minus Pass B)
    - B_not_C.bed                  (Pass B minus Pass C)
    - A_not_B_maf_ge_0.01.bed      (A-B subset with MAF >= 0.01)
    - A_not_B_maf_ge_0.05.bed      (A-B subset with MAF >= 0.05)
    - B_not_C_maf_ge_0.01.bed      (B-C subset with MAF >= 0.01)
    - B_not_C_maf_ge_0.05.bed      (B-C subset with MAF >= 0.05)

  Tables:
    - stage_summary.tsv            (Pass A/B/C site counts and attrition)
    - reason_summary.tsv           (per-reason drop counts)
    - debug_site_classification.tsv (per-site full table)
    - maf_summary.tsv              (MAF distribution summary from Pass A)
    - maf_threshold_counts.tsv     (retained sites at each MAF cutoff)
    - maf_bin_counts.tsv           (MAF histogram bins)

  Plots:
    - maf_histogram_0_0.02.png     (MAF zoom 0-0.02)
    - maf_histogram_0_0.2.png      (MAF zoom 0-0.2)
    - maf_histogram_log.png        (log-y full range)
    - maf_cumulative_retained.png  (cumulative retained vs MAF threshold)

Usage:
    python3 05B_classify_and_plot.py debug_meta.tsv

Exact vs approximate:
    - Classification: EXACT (set membership across ANGSD passes)
    - MAF/LRT values: EXACT (from ANGSD .mafs output)
    - nInd: EXACT if present in .mafs, otherwise >= minInd
    - FAIL_upstream count: EXACT (callable BED minus Pass A positions)
    - FAIL_not_callable: EXACT (complement of callable BED)
"""

import sys
import os
from collections import defaultdict

# ============================================================================
# Try to import matplotlib; if unavailable, skip plots
# ============================================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("[WARN] matplotlib not available; skipping plots. Tables + BEDs will still be produced.")


def load_meta(path):
    meta = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("key"):
                continue
            parts = line.split("\t", 1)
            if len(parts) == 2:
                meta[parts[0]] = parts[1]
    return meta


def load_mafs(path):
    """Load .mafs TSV into dict keyed by (chr, pos)."""
    sites = {}
    if not os.path.isfile(path):
        print(f"[WARN] Missing: {path}")
        return sites
    with open(path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < len(header):
                continue
            row = dict(zip(header, fields))
            chrom = row.get("chromo", "")
            pos = int(row.get("position", 0))
            sites[(chrom, pos)] = row
    return sites


def load_callable_intervals(bed_path, chrom_filter=None):
    intervals = defaultdict(list)
    if not os.path.isfile(bed_path):
        return intervals
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            c = parts[0]
            if chrom_filter and c != chrom_filter:
                continue
            intervals[c].append((int(parts[1]), int(parts[2])))
    for c in intervals:
        intervals[c].sort()
    return intervals


def callable_total_bp(intervals, chrom):
    return sum(e - s for s, e in intervals.get(chrom, []))


def get_float(row, *keys):
    for k in keys:
        if k in row:
            try:
                return float(row[k])
            except (ValueError, TypeError):
                pass
    return None


# ============================================================================
# MAIN
# ============================================================================
def main():
    if len(sys.argv) < 2:
        print("Usage: python3 05B_classify_and_plot.py debug_meta.tsv")
        sys.exit(1)

    meta = load_meta(sys.argv[1])
    chrom = meta["chrom"]
    chr_len = int(meta.get("chr_len", 0))
    start = int(meta.get("start", 0))
    end = int(meta.get("end", chr_len))
    snp_pval = meta["snp_pval"]
    maf_thresh = float(meta["maf"])
    outdir = meta["outdir"]
    tag = meta["tag"]
    minInd_str = meta.get("minInd", "?")

    print(f"[INFO] Tag: {tag}")
    print(f"[INFO] Region: {chrom}:{start}-{end}")

    # Load passes
    print("[INFO] Loading Pass A...")
    passA = load_mafs(meta["passA_mafs"])
    print(f"[INFO]   Pass A: {len(passA):,}")

    print("[INFO] Loading Pass B...")
    passB = load_mafs(meta["passB_mafs"])
    print(f"[INFO]   Pass B: {len(passB):,}")

    print("[INFO] Loading Pass C...")
    passC = load_mafs(meta["passC_mafs"])
    print(f"[INFO]   Pass C: {len(passC):,}")

    # Load callable
    print("[INFO] Loading callable BED...")
    callable_intervals = load_callable_intervals(meta["callable_bed"], chrom)
    callable_bp = callable_total_bp(callable_intervals, chrom)
    print(f"[INFO]   Callable bp: {callable_bp:,}")

    # ========================================================================
    # CLASSIFY Pass A sites
    # ========================================================================
    classified = {}
    for key, row in passA.items():
        if key in passC:
            classified[key] = "PASS_all_filters"
        elif key in passB:
            classified[key] = "FAIL_minMaf_only"
        else:
            maf_val = get_float(row, "unknownEM", "knownEM")
            if maf_val is not None and maf_val >= maf_thresh:
                classified[key] = "FAIL_snpPval_only"
            else:
                classified[key] = "FAIL_multiple_filters"

    cat_counts = defaultdict(int)
    for cat in classified.values():
        cat_counts[cat] += 1

    # ========================================================================
    # FIND UPSTREAM FAILURES (callable but not in Pass A)
    # ========================================================================
    print("[INFO] Finding upstream failures...")
    passA_positions = {p for (c, p) in passA if c == chrom}

    upstream_intervals = []
    upstream_count = 0
    run_start = None
    run_end = None

    for (s0, e0) in callable_intervals.get(chrom, []):
        for pos0 in range(s0, e0):
            pos1 = pos0 + 1
            if pos1 not in passA_positions:
                upstream_count += 1
                if run_start is None:
                    run_start = pos0
                    run_end = pos0 + 1
                elif pos0 == run_end:
                    run_end = pos0 + 1
                else:
                    upstream_intervals.append((run_start, run_end))
                    run_start = pos0
                    run_end = pos0 + 1
            else:
                if run_start is not None:
                    upstream_intervals.append((run_start, run_end))
                    run_start = None
    if run_start is not None:
        upstream_intervals.append((run_start, run_end))

    print(f"[INFO]   FAIL_upstream: {upstream_count:,} bp in {len(upstream_intervals):,} intervals")

    # Not-callable intervals
    not_callable = []
    if chr_len > 0:
        region_start = start if start > 0 else 0
        region_end = end if end > 0 else chr_len
        prev = region_start
        for (s0, e0) in callable_intervals.get(chrom, []):
            cs = max(s0, region_start)
            ce = min(e0, region_end)
            if cs >= ce:
                continue
            if cs > prev:
                not_callable.append((prev, cs))
            prev = ce
        if prev < region_end:
            not_callable.append((prev, region_end))
    not_callable_bp = sum(e - s for s, e in not_callable)

    # ========================================================================
    # WRITE BED TRACKS
    # ========================================================================
    print("[INFO] Writing BED tracks...")

    BED_COLORS = {
        "PASS_all_filters":      "0,150,0",
        "FAIL_minMaf_only":      "255,165,0",
        "FAIL_snpPval_only":     "255,0,0",
        "FAIL_multiple_filters": "139,0,139",
        "FAIL_upstream":         "100,100,100",
        "FAIL_not_callable":     "0,0,0",
    }

    def write_cat_bed(category, filename):
        path = os.path.join(outdir, filename)
        sites = sorted([(c, p) for (c, p), cat in classified.items() if cat == category], key=lambda x: x[1])
        with open(path, "w") as f:
            f.write(f'track name="{category}" color={BED_COLORS.get(category, "128,128,128")}\n')
            for (c, p) in sites:
                f.write(f"{c}\t{p-1}\t{p}\t{category}\n")
        print(f"[INFO]   {filename}: {len(sites):,}")

    write_cat_bed("PASS_all_filters", "PASS_all_filters.bed")
    write_cat_bed("FAIL_minMaf_only", "FAIL_minMaf_only.bed")
    write_cat_bed("FAIL_snpPval_only", "FAIL_snpPval_only.bed")
    write_cat_bed("FAIL_multiple_filters", "FAIL_multiple_filters.bed")

    # Interval BEDs
    def write_interval_bed(intervals_list, filename, name, color):
        path = os.path.join(outdir, filename)
        with open(path, "w") as f:
            f.write(f'track name="{name}" color={color}\n')
            for (s, e) in intervals_list:
                f.write(f"{chrom}\t{s}\t{e}\t{name}\n")
        print(f"[INFO]   {filename}: {len(intervals_list):,} intervals")

    write_interval_bed(upstream_intervals, "FAIL_upstream.bed", "FAIL_upstream", "100,100,100")
    write_interval_bed(not_callable, "FAIL_not_callable.bed", "FAIL_not_callable", "0,0,0")

    # Set-difference BEDs (Part 7)
    def write_setdiff_bed(setA, setB, filename, name, color, maf_cutoff=None):
        path = os.path.join(outdir, filename)
        diff_keys = []
        for key in setA:
            if key not in setB:
                if maf_cutoff is not None:
                    maf_val = get_float(setA[key], "unknownEM", "knownEM")
                    if maf_val is None or maf_val < maf_cutoff:
                        continue
                diff_keys.append(key)
        diff_keys.sort(key=lambda x: x[1])
        with open(path, "w") as f:
            f.write(f'track name="{name}" color={color}\n')
            for (c, p) in diff_keys:
                row = setA[(c, p)]
                maj = row.get("major", ".")
                mi = row.get("minor", ".")
                maf_v = get_float(row, "unknownEM", "knownEM")
                label = f"{maj}>{mi}_maf={maf_v:.4f}" if maf_v else f"{maj}>{mi}"
                f.write(f"{c}\t{p-1}\t{p}\t{label}\n")
        print(f"[INFO]   {filename}: {len(diff_keys):,}")

    write_setdiff_bed(passA, passB, "A_not_B.bed", "A_not_B", "200,50,50")
    write_setdiff_bed(passB, passC, "B_not_C.bed", "B_not_C", "200,150,50")
    write_setdiff_bed(passA, passB, "A_not_B_maf_ge_0.01.bed", "A_not_B_maf>=0.01", "200,80,80", 0.01)
    write_setdiff_bed(passA, passB, "A_not_B_maf_ge_0.05.bed", "A_not_B_maf>=0.05", "200,100,100", 0.05)
    write_setdiff_bed(passB, passC, "B_not_C_maf_ge_0.01.bed", "B_not_C_maf>=0.01", "200,180,80", 0.01)
    write_setdiff_bed(passB, passC, "B_not_C_maf_ge_0.05.bed", "B_not_C_maf>=0.05", "200,200,100", 0.05)

    # ========================================================================
    # WRITE TABLES (Part 7)
    # ========================================================================
    print("[INFO] Writing tables...")

    nA = len(passA)
    nB = len(passB)
    nC = len(passC)

    # Stage summary
    with open(os.path.join(outdir, "stage_summary.tsv"), "w") as f:
        f.write("stage\tn_sites\tdropped_vs_previous\tdropped_vs_initial\tfraction_initial\n")
        f.write(f"passA_permissive\t{nA}\t0\t0\t1.0000\n")
        f.write(f"passB_snpPval\t{nB}\t{nA-nB}\t{nA-nB}\t{nB/nA:.4f}\n") if nA > 0 else None
        f.write(f"passC_production\t{nC}\t{nB-nC}\t{nA-nC}\t{nC/nA:.4f}\n") if nA > 0 else None

    # Reason summary
    with open(os.path.join(outdir, "reason_summary.tsv"), "w") as f:
        f.write("reason\tn_sites_dropped\tfraction_initial\tfraction_callable\n")
        total_init = callable_bp if callable_bp > 0 else 1
        f.write(f"FAIL_not_callable\t{not_callable_bp}\tNA\t{not_callable_bp/total_init:.6f}\n")
        f.write(f"FAIL_upstream\t{upstream_count}\t{upstream_count/total_init:.6f}\tNA\n")
        for cat in ["FAIL_snpPval_only", "FAIL_minMaf_only", "FAIL_multiple_filters"]:
            n = cat_counts[cat]
            fi = n / nA if nA > 0 else 0
            f.write(f"{cat}\t{n}\t{fi:.6f}\t{n/total_init:.6f}\n")
        f.write(f"PASS_all_filters\t{cat_counts['PASS_all_filters']}\t{cat_counts['PASS_all_filters']/nA:.6f}\t{cat_counts['PASS_all_filters']/total_init:.6f}\n") if nA > 0 else None

    # Per-site classification TSV
    tsv_path = os.path.join(outdir, "debug_site_classification.tsv")
    with open(tsv_path, "w") as f:
        f.write("chr\tpos\tmajor\tminor\tnInd\tknownEM\tlrt\tpK_EM\tfail_reason\n")
        for key in sorted(classified.keys(), key=lambda x: x[1]):
            (c, p) = key
            cat = classified[key]
            row = passA[key]
            maj = row.get("major", ".")
            mi = row.get("minor", ".")
            nind = row.get("nInd", f">={minInd_str}")
            maf_v = get_float(row, "knownEM", "unknownEM")
            lrt_v = get_float(row, "LRT")  # may not be present in all .mafs
            pk = row.get("pK-EM", "NA")
            maf_s = f"{maf_v:.6f}" if maf_v is not None else "NA"
            lrt_s = f"{lrt_v:.4f}" if lrt_v is not None else "NA"
            f.write(f"{c}\t{p}\t{maj}\t{mi}\t{nind}\t{maf_s}\t{lrt_s}\t{pk}\t{cat}\n")
    print(f"[INFO]   debug_site_classification.tsv: {len(classified):,} rows")

    # ========================================================================
    # MAF ANALYSIS (Part 8) — from Pass A .mafs data
    # ========================================================================
    print("[INFO] Computing MAF distributions from Pass A...")

    mafs_all = []
    for key, row in passA.items():
        v = get_float(row, "knownEM", "unknownEM")
        if v is not None:
            mafs_all.append(v)

    mafs_all.sort()
    n_total = len(mafs_all)

    # MAF summary
    with open(os.path.join(outdir, "maf_summary.tsv"), "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"n_sites_passA\t{n_total}\n")
        if n_total > 0:
            f.write(f"maf_min\t{mafs_all[0]:.6f}\n")
            f.write(f"maf_max\t{mafs_all[-1]:.6f}\n")
            f.write(f"maf_median\t{mafs_all[n_total//2]:.6f}\n")
            f.write(f"maf_mean\t{sum(mafs_all)/n_total:.6f}\n")
            n_mono = sum(1 for m in mafs_all if m < 0.001)
            f.write(f"n_maf_below_0.001\t{n_mono}\n")
            f.write(f"n_maf_below_0.01\t{sum(1 for m in mafs_all if m < 0.01)}\n")
            f.write(f"n_maf_below_0.05\t{sum(1 for m in mafs_all if m < 0.05)}\n")

    # MAF threshold counts (cumulative retained)
    thresholds = [0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
    with open(os.path.join(outdir, "maf_threshold_counts.tsv"), "w") as f:
        f.write("maf_threshold\tn_retained\tfraction_passA\n")
        for t in thresholds:
            n = sum(1 for m in mafs_all if m >= t)
            f.write(f"{t:.3f}\t{n}\t{n/n_total:.6f}\n") if n_total > 0 else f.write(f"{t:.3f}\t0\t0\n")

    # MAF bin counts (fine bins)
    bin_width = 0.001
    bins = defaultdict(int)
    for m in mafs_all:
        b = int(m / bin_width) * bin_width
        bins[round(b, 4)] += 1

    with open(os.path.join(outdir, "maf_bin_counts.tsv"), "w") as f:
        f.write("maf_bin_start\tmaf_bin_end\tcount\n")
        for b in sorted(bins.keys()):
            f.write(f"{b:.4f}\t{b+bin_width:.4f}\t{bins[b]}\n")

    # ========================================================================
    # PLOTS (Part 8)
    # ========================================================================
    if HAS_MPL and n_total > 0:
        print("[INFO] Generating plots...")

        # Plot 1: MAF histogram 0-0.02 (zoom into rare variants)
        fig, ax = plt.subplots(figsize=(10, 5))
        maf_rare = [m for m in mafs_all if m <= 0.02]
        ax.hist(maf_rare, bins=200, color='steelblue', edgecolor='none', alpha=0.8)
        ax.set_xlabel("Minor Allele Frequency (MAF)")
        ax.set_ylabel("Number of sites")
        ax.set_title(f"MAF distribution (0–0.02) — Pass A, {tag}\nn={len(maf_rare):,} sites with MAF ≤ 0.02 out of {n_total:,}")
        ax.set_xlim(0, 0.02)
        ax.axvline(x=float(meta["maf"]), color='red', linestyle='--', label=f'minMaf={meta["maf"]}')
        ax.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "maf_histogram_0_0.02.png"), dpi=200)
        plt.close()

        # Plot 2: MAF histogram 0-0.2
        fig, ax = plt.subplots(figsize=(10, 5))
        maf_mid = [m for m in mafs_all if m <= 0.2]
        ax.hist(maf_mid, bins=200, color='steelblue', edgecolor='none', alpha=0.8)
        ax.set_xlabel("Minor Allele Frequency (MAF)")
        ax.set_ylabel("Number of sites")
        ax.set_title(f"MAF distribution (0–0.2) — Pass A, {tag}\nn={len(maf_mid):,} sites")
        ax.set_xlim(0, 0.2)
        ax.axvline(x=float(meta["maf"]), color='red', linestyle='--', label=f'minMaf={meta["maf"]}')
        ax.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "maf_histogram_0_0.2.png"), dpi=200)
        plt.close()

        # Plot 3: Full range log-y
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.hist(mafs_all, bins=500, color='steelblue', edgecolor='none', alpha=0.8)
        ax.set_yscale('log')
        ax.set_xlabel("Minor Allele Frequency (MAF)")
        ax.set_ylabel("Number of sites (log scale)")
        ax.set_title(f"MAF distribution (full range, log-y) — Pass A, {tag}\nn={n_total:,}")
        ax.set_xlim(0, 0.5)
        ax.axvline(x=float(meta["maf"]), color='red', linestyle='--', label=f'minMaf={meta["maf"]}')
        ax.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "maf_histogram_log.png"), dpi=200)
        plt.close()

        # Plot 4: Cumulative retained sites across MAF thresholds
        fig, ax = plt.subplots(figsize=(10, 5))
        thresholds_fine = [i * 0.001 for i in range(501)]
        retained = []
        for t in thresholds_fine:
            retained.append(sum(1 for m in mafs_all if m >= t))
        ax.plot(thresholds_fine, retained, color='steelblue', linewidth=1.5, label=f'SNP_pval={snp_pval}')
        ax.set_xlabel("MAF threshold")
        ax.set_ylabel("Number of retained sites")
        ax.set_title(f"Cumulative retained sites vs MAF threshold — Pass A, {tag}")
        ax.axvline(x=float(meta["maf"]), color='red', linestyle='--', alpha=0.7, label=f'minMaf={meta["maf"]}')
        ax.legend()
        ax.set_xlim(0, 0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "maf_cumulative_retained.png"), dpi=200)
        plt.close()

        print("[INFO] Plots saved.")
    elif not HAS_MPL:
        print("[WARN] Skipping plots (matplotlib not available)")

    # ========================================================================
    # PRINT SUMMARY
    # ========================================================================
    print("")
    print("=" * 72)
    print(f"  SUMMARY: {tag}  (SNP_pval={snp_pval}, MAF={maf_thresh})")
    print("=" * 72)
    if chr_len > 0:
        print(f"  Region:               {chrom}:{start}-{end}")
        print(f"  Not callable:         {not_callable_bp:>12,} bp")
        print(f"  Callable:             {callable_bp:>12,} bp")
    print(f"  FAIL_upstream:        {upstream_count:>12,} positions")
    print(f"  Pass A (eligible):    {nA:>12,}")
    print(f"    FAIL_snpPval_only:  {cat_counts['FAIL_snpPval_only']:>12,}")
    print(f"    FAIL_minMaf_only:   {cat_counts['FAIL_minMaf_only']:>12,}")
    print(f"    FAIL_multiple:      {cat_counts['FAIL_multiple_filters']:>12,}")
    print(f"    PASS_all_filters:   {cat_counts['PASS_all_filters']:>12,}")
    print(f"  Pass B (pval only):   {nB:>12,}")
    print(f"  Pass C (production):  {nC:>12,}")
    print("=" * 72)
    print(f"\n  Output: {outdir}/")
    print("[DONE]")


if __name__ == "__main__":
    main()
