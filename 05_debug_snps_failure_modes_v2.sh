###
# 05_debug_snps_failure_modes_v2.sh does this:
# Runs up to three progressive ANGSD passes on a region/chromosome to decompose
# SNP filter failure modes.  Supports restart-safe logic (skip existing outputs),
# region subsetting (--start/--end), selective pass execution (--run-only),
# and the new ANGSD -doFilterAudit patch for native audit counters.
#
# Pass A: permissive (SNP_pval=1, minMaf=0) — all upstream-eligible sites
# Pass B: SNP_pval only (no minMaf)           — sites passing LRT
# Pass C: production (SNP_pval + minMaf)       — final panel
#
# Usage:
#   bash 05_debug_snps_failure_modes_v2.sh --chrom C_gar_LG28
#   bash 05_debug_snps_failure_modes_v2.sh --chrom C_gar_LG28 --start 10108000 --end 10208000
#   bash 05_debug_snps_failure_modes_v2.sh --chrom C_gar_LG28 --snp-pval 1e-2
#   bash 05_debug_snps_failure_modes_v2.sh --chrom C_gar_LG28 --run-only passA
#   bash 05_debug_snps_failure_modes_v2.sh --chrom C_gar_LG28 --run-only classify
#   bash 05_debug_snps_failure_modes_v2.sh --chrom C_gar_LG28 --force
#
# For 100 kb debug-pack mode:
#   bash 05_debug_snps_failure_modes_v2.sh \
#     --chrom C_gar_LG28 --start 10108000 --end 10208000 \
#     --bamlist /path/to/debug_pack/meta/bamlist.txt \
#     --angsd /path/to/patched/angsd
###
#!/usr/bin/env bash
set -euo pipefail

# ========================= DEFAULT CONFIGURABLE PATHS ========================
BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
BAMLIST="${BASE}/bamlist.pp.samechr.tlenP99.filtered.txt"
PEST_SFS="${BASE}/popstruct_global/global_sfs/catfish.global.folded.sfs.mean.pest"
CALLABLE_SITES="${BASE}/fClaHyb_Gar_LG.mask_regions.normalACGT.renamed.1-based.angsd"
CALLABLE_BED="${BASE}/fClaHyb_Gar_LG.mask_regions.normalACGT.renamed.0based.bed"
FAI="${REF}.fai"

# ANGSD binary — set to patched version if available
ANGSD_BIN="angsd"

# ========================= DEFAULT FILTER SETTINGS ===========================
MINQ=25
MINMAPQ=25
MINDEPTHIND=3
MAXDEPTHIND=57
MININD=200
SNP_PVAL="1e-6"
MAF="0.05"
REMOVE_FLAGS=1
UNIQUE_ONLY=1
ONLY_PROPER_PAIRS=1
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ========================= PARSE ARGUMENTS ===================================
CHROM=""
RF_FILE=""
START=""
END=""
RUN_ONLY="all"
FORCE=0
DO_AUDIT=1   # 0=no audit, 1=summary, 2=per-site TSV

usage() {
    cat <<EOF
Usage:
  $0 --chrom NAME [options]
  $0 --rf FILE [options]

Options:
  --chrom NAME       Chromosome name
  --rf FILE          Path to .rf.txt
  --start POS        Start position (1-based, for region subset)
  --end POS          End position (1-based, for region subset)
  --snp-pval VAL     SNP p-value threshold (default: 1e-6)
  --maf VAL          MAF threshold (default: 0.05)
  --minind VAL       minInd filter (default: 200)
  --threads N        Threads (default: SLURM_CPUS_PER_TASK or 8)
  --bamlist FILE     BAM list (default: production bamlist)
  --callable FILE    ANGSD callable sites file
  --angsd PATH       Path to ANGSD binary
  --run-only STAGE   Run only: passA|passB|passC|classify|all (default: all)
  --force            Force rerun even if outputs exist
  --audit LEVEL      0=off, 1=summary, 2=per-site TSV (default: 1)
  -h|--help          Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrom)      CHROM="$2"; shift 2 ;;
        --rf)         RF_FILE="$2"; shift 2 ;;
        --start)      START="$2"; shift 2 ;;
        --end)        END="$2"; shift 2 ;;
        --snp-pval)   SNP_PVAL="$2"; shift 2 ;;
        --maf)        MAF="$2"; shift 2 ;;
        --minind)     MININD="$2"; shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        --bamlist)    BAMLIST="$2"; shift 2 ;;
        --callable)   CALLABLE_SITES="$2"; shift 2 ;;
        --angsd)      ANGSD_BIN="$2"; shift 2 ;;
        --run-only)   RUN_ONLY="$2"; shift 2 ;;
        --force)      FORCE=1; shift ;;
        --audit)      DO_AUDIT="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

[[ -n "$CHROM" || -n "$RF_FILE" ]] || { echo "[ERROR] Need --chrom or --rf" >&2; usage >&2; exit 1; }

# ========================= DERIVE TAG AND PATHS ==============================
if [[ -n "$CHROM" ]]; then
    TAG="$CHROM"
else
    TAG="$(basename "$RF_FILE" .rf.txt)"
    CHROM="$(head -1 "$RF_FILE")"
fi

# Region tag for output directory
if [[ -n "$START" && -n "$END" ]]; then
    REGION_TAG="${TAG}_${START}_${END}"
else
    REGION_TAG="$TAG"
fi

PVAL_TAG="${SNP_PVAL}"
MAF_TAG="${MAF}"
OUTDIR="${BASE}/debug_snp_failure/debug_${REGION_TAG}_pval${PVAL_TAG}_maf${MAF_TAG}"
mkdir -p "${OUTDIR}"/{passA,passB,passC,classify,logs}

# Create RF file
if [[ -z "$RF_FILE" ]]; then
    if [[ -n "$START" && -n "$END" ]]; then
        RF_FILE="${OUTDIR}/${REGION_TAG}.rf.txt"
        echo "${CHROM}:${START}-${END}" > "$RF_FILE"
    else
        RF_FILE="${OUTDIR}/${TAG}.rf.txt"
        echo "$CHROM" > "$RF_FILE"
    fi
fi

CHR_LEN=$(awk -v c="$CHROM" '$1==c {print $2}' "$FAI")

echo "========================================================================"
echo " SNP failure-mode debugging v2"
echo "========================================================================"
echo "[INFO] REGION_TAG    = $REGION_TAG"
echo "[INFO] CHROM         = $CHROM"
[[ -n "$START" ]] && echo "[INFO] START         = $START"
[[ -n "$END" ]]   && echo "[INFO] END           = $END"
echo "[INFO] RF_FILE       = $RF_FILE"
echo "[INFO] SNP_PVAL      = $SNP_PVAL"
echo "[INFO] MAF           = $MAF"
echo "[INFO] MININD        = $MININD"
echo "[INFO] OUTDIR        = $OUTDIR"
echo "[INFO] ANGSD         = $ANGSD_BIN"
echo "[INFO] RUN_ONLY      = $RUN_ONLY"
echo "[INFO] FORCE         = $FORCE"
echo "[INFO] DO_AUDIT      = $DO_AUDIT"
echo ""

# ========================= HELPER FUNCTIONS ==================================

is_valid_mafs() {
    local f="$1"
    [[ -f "$f" ]] && [[ -s "$f" ]] && zcat "$f" | head -1 | grep -q "^chromo"
}

should_run_pass() {
    local pass="$1"
    local mafs_file="$2"
    if [[ "$RUN_ONLY" != "all" && "$RUN_ONLY" != "$pass" ]]; then
        echo "[SKIP] $pass (--run-only=$RUN_ONLY)"
        return 1
    fi
    if [[ "$FORCE" -eq 0 ]] && is_valid_mafs "$mafs_file"; then
        echo "[SKIP] $pass — output exists and valid: $mafs_file"
        return 1
    fi
    return 0
}

# ========================= COMMON ANGSD FLAGS ================================
COMMON_FLAGS=(
    -b "$BAMLIST" -ref "$REF"
    -GL 1 -doMajorMinor 1 -doMaf 1
    -pest "$PEST_SFS"
    -doCounts 1 -doDepth 1
    -minQ "$MINQ" -minMapQ "$MINMAPQ" -baq 1 -C 50
    -setMinDepthInd "$MINDEPTHIND" -setMaxDepthInd "$MAXDEPTHIND"
    -minInd "$MININD"
    -remove_bads "$REMOVE_FLAGS"
    -uniqueOnly "$UNIQUE_ONLY"
    -only_proper_pairs "$ONLY_PROPER_PAIRS"
    -rf "$RF_FILE"
    -sites "$CALLABLE_SITES"
    -P "$THREADS"
)

# Add audit flags if patched ANGSD supports them
AUDIT_FLAGS=()
if [[ "$DO_AUDIT" -ge 1 ]]; then
    AUDIT_FLAGS=(-doFilterAudit "$DO_AUDIT")
fi

# =============================================================================
# PASS A: MAXIMALLY PERMISSIVE
# =============================================================================
PASS_A_PREFIX="${OUTDIR}/passA/catfish.${REGION_TAG}.passA"
PASS_A_MAFS="${PASS_A_PREFIX}.mafs.gz"

if should_run_pass "passA" "$PASS_A_MAFS"; then
    echo "========================================================================"
    echo " PASS A: Permissive (SNP_pval=1, minMaf=0)"
    echo "========================================================================"

    AUDIT_A_FLAGS=("${AUDIT_FLAGS[@]}")
    if [[ "$DO_AUDIT" -ge 2 ]]; then
        AUDIT_A_FLAGS+=(-filterAuditFile "${OUTDIR}/passA/audit_passA.tsv")
    fi

    ${ANGSD_BIN} "${COMMON_FLAGS[@]}" \
        -SNP_pval 1 -minMaf 0 \
        "${AUDIT_A_FLAGS[@]}" \
        -out "$PASS_A_PREFIX" \
        > "${OUTDIR}/logs/passA.log" 2>&1 || true

    if is_valid_mafs "$PASS_A_MAFS"; then
        echo "[INFO] Pass A done: $(zcat "$PASS_A_MAFS" | tail -n+2 | wc -l) sites"
    else
        echo "[ERROR] Pass A failed — check ${OUTDIR}/logs/passA.log"
        # If audit flags caused failure, retry without them
        echo "[INFO] Retrying Pass A without audit flags..."
        ${ANGSD_BIN} "${COMMON_FLAGS[@]}" \
            -SNP_pval 1 -minMaf 0 \
            -out "$PASS_A_PREFIX" \
            > "${OUTDIR}/logs/passA.log" 2>&1
        echo "[INFO] Pass A (no audit): $(zcat "$PASS_A_MAFS" | tail -n+2 | wc -l) sites"
    fi
fi

# =============================================================================
# PASS B: SNP_PVAL ONLY
# =============================================================================
PASS_B_PREFIX="${OUTDIR}/passB/catfish.${REGION_TAG}.passB"
PASS_B_MAFS="${PASS_B_PREFIX}.mafs.gz"

if should_run_pass "passB" "$PASS_B_MAFS"; then
    echo "========================================================================"
    echo " PASS B: SNP_pval=${SNP_PVAL} only (no minMaf)"
    echo "========================================================================"

    AUDIT_B_FLAGS=("${AUDIT_FLAGS[@]}")
    if [[ "$DO_AUDIT" -ge 2 ]]; then
        AUDIT_B_FLAGS+=(-filterAuditFile "${OUTDIR}/passB/audit_passB.tsv")
    fi

    ${ANGSD_BIN} "${COMMON_FLAGS[@]}" \
        -SNP_pval "$SNP_PVAL" -minMaf 0 \
        "${AUDIT_B_FLAGS[@]}" \
        -out "$PASS_B_PREFIX" \
        > "${OUTDIR}/logs/passB.log" 2>&1 || true

    if is_valid_mafs "$PASS_B_MAFS"; then
        echo "[INFO] Pass B done: $(zcat "$PASS_B_MAFS" | tail -n+2 | wc -l) sites"
    else
        echo "[ERROR] Pass B failed — check ${OUTDIR}/logs/passB.log"
        ${ANGSD_BIN} "${COMMON_FLAGS[@]}" \
            -SNP_pval "$SNP_PVAL" -minMaf 0 \
            -out "$PASS_B_PREFIX" \
            > "${OUTDIR}/logs/passB.log" 2>&1
        echo "[INFO] Pass B (no audit): $(zcat "$PASS_B_MAFS" | tail -n+2 | wc -l) sites"
    fi
fi

# =============================================================================
# PASS C: PRODUCTION
# =============================================================================
PASS_C_PREFIX="${OUTDIR}/passC/catfish.${REGION_TAG}.passC"
PASS_C_MAFS="${PASS_C_PREFIX}.mafs.gz"

if should_run_pass "passC" "$PASS_C_MAFS"; then
    echo "========================================================================"
    echo " PASS C: Production (SNP_pval=${SNP_PVAL}, minMaf=${MAF})"
    echo "========================================================================"

    AUDIT_C_FLAGS=("${AUDIT_FLAGS[@]}")
    if [[ "$DO_AUDIT" -ge 2 ]]; then
        AUDIT_C_FLAGS+=(-filterAuditFile "${OUTDIR}/passC/audit_passC.tsv")
    fi

    ${ANGSD_BIN} "${COMMON_FLAGS[@]}" \
        -SNP_pval "$SNP_PVAL" -minMaf "$MAF" \
        "${AUDIT_C_FLAGS[@]}" \
        -out "$PASS_C_PREFIX" \
        > "${OUTDIR}/logs/passC.log" 2>&1 || true

    if is_valid_mafs "$PASS_C_MAFS"; then
        echo "[INFO] Pass C done: $(zcat "$PASS_C_MAFS" | tail -n+2 | wc -l) sites"
    else
        echo "[ERROR] Pass C failed — check ${OUTDIR}/logs/passC.log"
        ${ANGSD_BIN} "${COMMON_FLAGS[@]}" \
            -SNP_pval "$SNP_PVAL" -minMaf "$MAF" \
            -out "$PASS_C_PREFIX" \
            > "${OUTDIR}/logs/passC.log" 2>&1
        echo "[INFO] Pass C (no audit): $(zcat "$PASS_C_MAFS" | tail -n+2 | wc -l) sites"
    fi
fi

# =============================================================================
# PREPARE METADATA FOR CLASSIFIER
# =============================================================================
if [[ "$RUN_ONLY" == "all" || "$RUN_ONLY" == "classify" ]]; then

    # Decompress for classifier
    for PASS in passA passB passC; do
        MAFS="${OUTDIR}/${PASS}/catfish.${REGION_TAG}.${PASS}.mafs.gz"
        TSV="${OUTDIR}/classify/${PASS}.mafs.tsv"
        if [[ -f "$MAFS" ]] && { [[ ! -f "$TSV" ]] || [[ "$FORCE" -eq 1 ]]; }; then
            zcat "$MAFS" > "$TSV"
            echo "[INFO] Decompressed ${PASS}: $(tail -n+2 "$TSV" | wc -l) sites"
        fi
    done

    # Extract callable BED for this chromosome
    CALLABLE_CHR_BED="${OUTDIR}/classify/${REGION_TAG}.callable.bed"
    if [[ ! -f "$CALLABLE_CHR_BED" ]] || [[ "$FORCE" -eq 1 ]]; then
        if [[ -n "$START" && -n "$END" ]]; then
            awk -v c="$CHROM" -v s="$((START-1))" -v e="$END" \
                '$1==c && $3>s && $2<e {a=$2; b=$3; if(a<s)a=s; if(b>e)b=e; print $1"\t"a"\t"b}' \
                "$CALLABLE_BED" > "$CALLABLE_CHR_BED"
        else
            grep -w "^${CHROM}" "$CALLABLE_BED" > "$CALLABLE_CHR_BED" || true
        fi
        echo "[INFO] Callable intervals: $(wc -l < "$CALLABLE_CHR_BED")"
    fi

    # Write metadata
    cat > "${OUTDIR}/classify/debug_meta.tsv" <<METAEOF
key	value
tag	${REGION_TAG}
chrom	${CHROM}
chr_len	${CHR_LEN}
start	${START:-0}
end	${END:-${CHR_LEN}}
snp_pval	${SNP_PVAL}
maf	${MAF}
minInd	${MININD}
callable_bed	${CALLABLE_CHR_BED}
passA_mafs	${OUTDIR}/classify/passA.mafs.tsv
passB_mafs	${OUTDIR}/classify/passB.mafs.tsv
passC_mafs	${OUTDIR}/classify/passC.mafs.tsv
outdir	${OUTDIR}/classify
METAEOF

    echo ""
    echo "========================================================================"
    echo " Run classifier:"
    echo "========================================================================"
    echo ""
    echo "  python3 05B_classify_and_plot.py ${OUTDIR}/classify/debug_meta.tsv"
    echo ""
fi

echo "[DONE] 05_debug_snps_failure_modes_v2.sh"
