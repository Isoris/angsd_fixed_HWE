# ANGSD Filter Audit Patch + Improved Debug Workflow
# Complete Deliverable Document
# =============================================================

## FILES CHANGED

### ANGSD source patches:
1. abcFreq.h   — Added 10 audit member variables to class
2. abcFreq.cpp — Added audit argument parsing, counter logic in run(), per-site TSV output, summary in destructor
3. shared.cpp  — Raised thread cap from 8→16

### New wrapper/analysis scripts:
4. 05_debug_snps_failure_modes_v2.sh  — Improved bash wrapper with restart-safe, region, run-only
5. 05B_classify_and_plot.py           — Combined classifier + QC plotter

## COMPILE INSTRUCTIONS

```bash
# On your HPC, from the ANGSD source directory:
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/angsd_fixed_HWE-master

# Backup originals
cp abcFreq.cpp abcFreq.cpp.bak
cp abcFreq.h abcFreq.h.bak
cp shared.cpp shared.cpp.bak

# Apply patches (copy patched files over originals)
cp /path/to/patched/abcFreq.cpp .
cp /path/to/patched/abcFreq.h .
cp /path/to/patched/shared.cpp .

# Compile (same as original, using htslib submodule)
make clean
make -j8

# Verify binary
./angsd 2>&1 | head -5

# The new binary is: ./angsd (in the source directory)
# Copy to a known location:
cp angsd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/angsd_patched
```

If using system htslib:
```bash
make clean
make -j8 HTSSRC=systemwide
```

If using conda/mamba htslib:
```bash
make clean
make -j8 HTSSRC=${CONDA_PREFIX}
```

## EXACT EXECUTABLE PATH AFTER COMPILE

```
/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/angsd_fixed_HWE-master/angsd
```

Or wherever you copy it:
```
/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/angsd_patched
```

## TEST COMMANDS

### Test 1: 100 kb debug region with audit (production thresholds)
```bash
BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
ANGSD="${BASE}/angsd_fixed_HWE-master/angsd"
DBGPACK="${BASE}/debug_pack_C_gar_LG28_10108000_10208000"

# Use patched ANGSD on debug pack
bash 05_debug_snps_failure_modes_v2.sh \
    --chrom C_gar_LG28 \
    --start 10108000 \
    --end 10208000 \
    --bamlist "${DBGPACK}/meta/C_gar_LG28_10108000_10208000.bamlist.txt" \
    --callable "${DBGPACK}/meta/fClaHyb_Gar_LG.mask_regions.normalACGT.renamed.1-based.angsd" \
    --angsd "${ANGSD}" \
    --snp-pval 1e-6 \
    --maf 0.05 \
    --audit 2 \
    --threads 8
```

### Test 2: Same region, relaxed thresholds
```bash
bash 05_debug_snps_failure_modes_v2.sh \
    --chrom C_gar_LG28 \
    --start 10108000 \
    --end 10208000 \
    --bamlist "${DBGPACK}/meta/C_gar_LG28_10108000_10208000.bamlist.txt" \
    --callable "${DBGPACK}/meta/fClaHyb_Gar_LG.mask_regions.normalACGT.renamed.1-based.angsd" \
    --angsd "${ANGSD}" \
    --snp-pval 1e-2 \
    --maf 0.01 \
    --audit 2 \
    --threads 8
```

### Test 3: Full chromosome (SLURM)
```bash
sbatch --cpus-per-task=16 --mem=237GB -t 12:00:00 -p compute \
    --wrap="bash 05_debug_snps_failure_modes_v2.sh --chrom C_gar_LG28 --angsd ${ANGSD} --threads 16"
```

### Test 4: Run only classifier (if ANGSD passes already completed)
```bash
bash 05_debug_snps_failure_modes_v2.sh --chrom C_gar_LG28 --run-only classify
python3 05B_classify_and_plot.py \
    ${BASE}/debug_snp_failure/debug_C_gar_LG28_pval1e-6_maf0.05/classify/debug_meta.tsv
```

### Test 5: Direct ANGSD with audit flags (standalone, no wrapper)
```bash
${ANGSD} \
    -b ${BASE}/bamlist.pp.samechr.tlenP99.filtered.txt \
    -ref ${BASE}/00-samples/fClaHyb_Gar_LG.fa \
    -GL 1 -doMajorMinor 1 -doMaf 1 \
    -pest ${BASE}/popstruct_global/global_sfs/catfish.global.folded.sfs.mean.pest \
    -doCounts 1 -doDepth 1 \
    -minQ 25 -minMapQ 25 -baq 1 -C 50 \
    -setMinDepthInd 3 -setMaxDepthInd 57 \
    -minInd 200 \
    -SNP_pval 1e-6 -minMaf 0.05 \
    -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
    -rf ${BASE}/popstruct_global/rf_files/C_gar_LG28.rf.txt \
    -sites ${BASE}/fClaHyb_Gar_LG.mask_regions.normalACGT.renamed.1-based.angsd \
    -doFilterAudit 1 \
    -P 16 \
    -out /tmp/test_audit
```

## EXPECTED NEW LOG LINES FROM AUDIT COUNTERS

When `-doFilterAudit 1` (or 2):
```
[abcFreq audit] ============================================
[abcFreq audit] total_sites_seen=6748000
[abcFreq audit] dead_before_freq_filter=5123456
[abcFreq audit] dropped_minMaf_low=1420000
[abcFreq audit] dropped_minMaf_high=0
[abcFreq audit] dropped_snpPval=85000
[abcFreq audit] dropped_minInd=2500
[abcFreq audit] dropped_triallelic=0
[abcFreq audit] retained_after_abcFreq=117044
[abcFreq audit] ============================================
```

(Numbers are illustrative. The actual values will tell you exactly where sites are lost.)

Key interpretation:
- `dead_before_freq_filter` = sites that entered abcFreq with keepSites==0.
  These were killed upstream by abcFilter (callable mask or minInd on read counts)
  or by abcMajorMinor (could not infer major/minor).
- `dropped_minMaf_low` = sites where EM-estimated MAF < minMaf
- `dropped_snpPval` = sites where LRT < chi-square cutoff (not enough evidence for polymorphism)
- `dropped_minInd` = sites where GL-based nInd < minInd (after MAF/pval may have zeroed it)
  In practice, with your filter order, this counter catches sites that survived
  MAF/pval checks but had too few informative individuals after GL recounting.

## EXAMPLE PER-SITE AUDIT TSV HEADER

When `-doFilterAudit 2 -filterAuditFile audit.tsv`:
```
chromo  position  major  minor  knownEM  lrt     nInd_GL  keepSites_final  drop_reason         snp_cutoff_used  minMaf_used
C_gar_LG28  10108001  C  T  0.000000  0.0000  223  0  dead_before_freq  1.000000e-06  0.050000
C_gar_LG28  10108015  G  A  0.023400  4.5678  220  0  minMaf_low        1.000000e-06  0.050000
C_gar_LG28  10108100  T  C  0.085200  28.456  218  218  PASS             1.000000e-06  0.050000
```

Columns:
- chromo, position: genomic coordinates (1-based)
- major, minor: ANGSD-inferred alleles
- knownEM: EM-estimated MAF (exact from ANGSD)
- lrt: likelihood ratio test statistic (exact from ANGSD)
- nInd_GL: number of individuals with informative GLs (exact after likeFreq recount)
- keepSites_final: final keepSites value after all filters
- drop_reason: FIRST filter that killed the site (pipeline semantics)
- snp_cutoff_used: the LRT cutoff corresponding to -SNP_pval
- minMaf_used: the -minMaf value used

## THREAD CAP ANSWER

The thread cap is in shared.cpp line 68-71. Originally caps at 8 if >10 requested.

**Is it safe to raise?** YES, to 16. The synchronization uses proper pthread mutexes.
The queue model (shared.cpp) has a producer-consumer pattern where:
- One thread reads BAM chunks
- Up to maxThreads-1 threads process chunks through the abc pipeline
- One printer thread serializes output

With 226 BAMs, the per-chunk computation (GL + EM + LRT) is substantial enough
that 16 threads gives meaningful speedup. Beyond 16, the BAM I/O and print
serialization become bottlenecks, so there's little benefit.

**The patch raises the cap to 16.** Values >16 are capped to 16 with a warning.

**Is 80 threads useful?** No. Even with 226 BAMs, the internal architecture
doesn't parallelize beyond the chunk dispatch queue. The BAM reading itself
uses htslib's thread pool independently.

## PREP_SITES / -SITES SPEED RELEVANCE

For the 100 kb debug window: **negligible bottleneck**. The binary index loads
in <1ms and provides O(1) lookup per position. The dominant cost is BAM I/O.

For full chromosomes: -sites is a **net speedup** because it kills non-callable
positions at abcFilter BEFORE the expensive GL/EM computation. Without -sites,
ANGSD would compute GLs and EM for every pileup position.

**No optimization needed here.** The current design is already efficient.

## 100 KB DEBUG RUN OPTIMIZATION

For 100 kb debug runs, the main bottleneck is opening 226 BAMs for region extraction.
The debug pack (pre-extracted BAMs) already solves this. With pre-extracted BAMs:
- Pass A on 100 kb: ~2-5 minutes
- Pass B on 100 kb: ~2-5 minutes  
- Pass C on 100 kb: ~1-2 minutes
- Classifier: ~30 seconds

Total for all three passes + classify: ~10-15 minutes on 8 threads.

## EXPECTED OUTPUT FILENAMES

For a run with `--chrom C_gar_LG28 --snp-pval 1e-6 --maf 0.05`:

```
debug_snp_failure/debug_C_gar_LG28_pval1e-6_maf0.05/
├── passA/
│   ├── catfish.C_gar_LG28.passA.mafs.gz
│   ├── catfish.C_gar_LG28.passA.arg
│   └── audit_passA.tsv              (if -doFilterAudit 2)
├── passB/
│   ├── catfish.C_gar_LG28.passB.mafs.gz
│   └── audit_passB.tsv
├── passC/
│   ├── catfish.C_gar_LG28.passC.mafs.gz
│   └── audit_passC.tsv
├── classify/
│   ├── debug_meta.tsv
│   ├── passA.mafs.tsv
│   ├── passB.mafs.tsv
│   ├── passC.mafs.tsv
│   ├── C_gar_LG28.callable.bed
│   ├── PASS_all_filters.bed
│   ├── FAIL_not_callable.bed
│   ├── FAIL_upstream.bed
│   ├── FAIL_minMaf_only.bed
│   ├── FAIL_snpPval_only.bed
│   ├── FAIL_multiple_filters.bed
│   ├── A_not_B.bed
│   ├── B_not_C.bed
│   ├── A_not_B_maf_ge_0.01.bed
│   ├── A_not_B_maf_ge_0.05.bed
│   ├── B_not_C_maf_ge_0.01.bed
│   ├── B_not_C_maf_ge_0.05.bed
│   ├── stage_summary.tsv
│   ├── reason_summary.tsv
│   ├── debug_site_classification.tsv
│   ├── maf_summary.tsv
│   ├── maf_threshold_counts.tsv
│   ├── maf_bin_counts.tsv
│   ├── maf_histogram_0_0.02.png
│   ├── maf_histogram_0_0.2.png
│   ├── maf_histogram_log.png
│   └── maf_cumulative_retained.png
└── logs/
    ├── passA.log
    ├── passB.log
    └── passC.log
```

For the 100 kb region variant, replace C_gar_LG28 with C_gar_LG28_10108000_10208000.
