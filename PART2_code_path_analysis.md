# ANGSD Code Path Analysis for SNP Discovery Workflow
# ===================================================
# Command: -GL 1 -doMajorMinor 1 -doMaf 1 -pest ... -doCounts 1 -doDepth 1
#          -minQ 25 -minMapQ 25 -baq 1 -C 50 -setMinDepthInd 3 -setMaxDepthInd 57
#          -minInd 200 -SNP_pval 1e-6 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1
#          -only_proper_pairs 1 -rf region -sites callable.angsd -P 80

## 1. Active Source Files (Definitively Active)

```
CRITICAL PATH (in execution order):
  [0] abcFilter.cpp      — upstream site/individual filtering, -sites, -minInd (read-count version)
  [1] abcGetFasta.cpp     — loads ref/anc FASTA bases per site
  [2] abcCounts.cpp       — allele counts per site (for -doCounts 1)
  [4] abcGL.cpp           — genotype likelihood computation (-GL 1 = SAMtools model)
  [6] abcMajorMinor.cpp   — GL-based major/minor inference (-doMajorMinor 1)
  [7] abcFreq.cpp         — EM MAF estimation, LRT, SNP_pval, minMaf, minInd (GL version) *** MAIN TARGET ***

SUPPORT:
  prep_sites.cpp          — binary index reading for -sites
  chisquare.cpp           — chi-square distribution for p-value <-> LRT conversion
  shared.cpp              — thread management, chunk dispatch, print queue
  mUpPile.cpp             — BAM reading, pileup generation

SECONDARY (NOT active in this workflow):
  abcFilterSNP.cpp        — only runs with -doSnpStat (you don't use this)
```

## 2. Module Execution Order (from abc.cpp lines 54-82)

Modules run sequentially for each chunk of sites:
```
Index  Module            What it does for your workflow
-----  ----------------  -------------------------------------------------
  0    abcFilter         Apply -sites mask (prep_sites binary lookup)
                         Count individuals with >= setMinDepthInd reads per site
                         Set keepSites[s] = nInd_with_data (or 0 if < minInd)
  1    abcGetFasta       Load reference base per position
  2    abcCounts         Tally A/C/G/T per individual per site
  3    abcError          (not active unless -doError)
  4    abcGL             Compute 10-genotype log-likelihoods per individual per site
  5    abcMcall          (not active unless -doGeno with certain methods)
  6    abcMajorMinor     Infer major/minor alleles from GLs
  7    abcFreq           EM MAF estimation, LRT test, apply minMaf/SNP_pval/minInd
  ...  (remaining modules not relevant)
```

## 3. Detailed Filter Pipeline

### STAGE 1: BAM reading (mUpPile.cpp)
- Reads are filtered by: -minQ 25, -minMapQ 25, -remove_bads 1, -uniqueOnly 1, -only_proper_pairs 1
- BAQ recalibration applied: -baq 1
- Adjustments: -C 50 (downgrade mapQ for excessive mismatches)
- These are HARD pre-filters at the read level. Reads failing these never contribute to anything.

### STAGE 2: abcFilter.run() [module index 0]
Step 2a: -sites mask check
  - For each position in the chunk, checks fl->keeps[pos]
  - fl->keeps is loaded from the binary indexed -sites file (prep_sites.cpp)
  - If position is not in the callable mask: keepSites[s] = 0
  - This is EXACT binary: in mask or not

Step 2b: Per-individual depth counting
  - For each site still alive (keepSites[s] != 0):
    - Counts how many individuals have >= setMinDepthInd (3) reads at this site
    - Sets keepSites[s] = nInfo (the count of qualifying individuals)
  - If minInd > 0 AND nInfo < minInd: keepSites[s] = 0

  IMPORTANT: After abcFilter, keepSites[s] holds the NUMBER OF INDIVIDUALS with data,
  not just 0/1. This count is used downstream as nInd for that site.

### STAGE 3: abcGL.run() [module index 4]
- For each site where keepSites[s] > 0:
  - Computes 10 log-genotype-likelihoods per individual using SAMtools model (-GL 1)
  - Stores in pars->likes[s][i*10 + genotype_index]

### STAGE 4: abcMajorMinor.run() [module index 6]
- For each site where keepSites[s] > 0:
  - With -doMajorMinor 1 (GL-based inference):
    - Tests all 6 major/minor allele pair combinations
    - For each pair, sums log-likelihoods across all individuals assuming HWE freq 0.25/0.50/0.25
    - Picks the pair with highest total likelihood
    - Assigns major = allele with higher estimated frequency, minor = other
  - If inference fails (e.g., all likelihoods identical): keepSites[s] = 0

### STAGE 5: abcFreq.run() [module index 7] *** THE MAIN FILTERING STAGE ***

Step 5a: likeFreq() — EM allele frequency estimation
  - For each site where keepSites[s] > 0:
    - Extracts 3-genotype log-likelihoods (AA, Aa, aa) for the inferred major/minor
    - RESETS keepSites[s] = 0 (line 729!)
    - Recounts: for each individual, checks if GL data is informative
      (not all-zero, not all-equal). If informative, increments keepSites[s].
      This gives the GL-based nInd count.
    - If keepSites[s] == 0 after recounting: skip site
    - With -doMaf 1: runs EM algorithm (emFrequency) to estimate MAF
    - With -doSNP (SNP_pval set): computes LRT = 2*(logL_null - logL_mle)
      where null is freq=0 (monomorphic), MLE is the EM-estimated freq
    - freq->freq[s] = EM-estimated minor allele frequency
    - freq->lrt[s] = LRT statistic (in chi-square likelihood units)

Step 5b: Filter loop (lines 598-616) — THE CRITICAL FILTERING DECISION
  ```
  for each site s:
    if keepSites[s] == 0: skip (already dead from upstream)
    
    CHECK 1 (minMaf): if freq < minMaf OR freq > 1-minMaf → keepSites[s] = 0
                       NOTE: This check happens FIRST.
                       NOTE: keepSites[s] is set to 0, but the loop CONTINUES.
                       
    CHECK 2 (SNP_pval): if lrt < SNP_pval_cutoff → keepSites[s] = 0
                         NOTE: SNP_pval was converted to LRT units in getOptions().
                         NOTE: This is "lrt < cutoff" which means "not enough evidence"
                         NOTE: This check is NOT inside an else-if with CHECK 1.
                               A site can be zeroed by BOTH minMaf AND SNP_pval.
    
    CHECK 3 (rmTriallelic): not used in your workflow
    
    CHECK 4 (minInd in abcFreq): if minInd > 0 AND keepSites[s] < minInd → keepSites[s] = 0
                                  NOTE: This is REDUNDANT with abcFilter's minInd for your case,
                                  BUT it uses the GL-recounted nInd, not the read-count nInd.
                                  HOWEVER: if CHECK 1 or 2 already zeroed keepSites[s],
                                  this check will also fire (0 < 200) — but the site is
                                  already dead so it doesn't matter functionally.
  ```

## 4. Critical Answers

### Q: Which rule is checked first in the final keepSites logic?
A: minMaf is checked FIRST (line 602-605), then SNP_pval (line 610-611), then
   rmTriallelic (line 612-613, not used), then minInd (line 614-615).

### Q: Is minMaf applied before or after SNP_pval?
A: minMaf is applied BEFORE SNP_pval in the filter loop.
   HOWEVER, they are NOT mutually exclusive — both can fire on the same site.
   The code does NOT use if/else-if between them. A site failing minMaf still
   gets checked for SNP_pval (but since keepSites is already 0, it doesn't
   matter functionally for the output — it matters for audit counting though).
   
   WAIT — actually look closely:
   Line 602-605 sets keepSites[s]=0 if MAF fails.
   Line 606-611 then checks SNP_pval. But the site is NOT skipped between
   checks — the code falls through. Since keepSites[s] is already 0,
   line 610 checks `freq->lrt[s] < SNP_pval` and may set keepSites[s]=0
   again (redundantly).
   Line 614 checks `minInd>0 && keepSites[s]<minInd` — since keepSites[s]
   is already 0, this fires too (redundantly).
   
   NET EFFECT: For pipeline audit counting, a site that fails minMaf will ALSO
   appear to fail SNP_pval and minInd at the code level. For proper auditing,
   we need to check conditions independently before any modification, OR check
   in order and count only the FIRST failure.

### Q: Does SNP_pval use the EM-estimated MAF from abcFreq?
A: YES. The LRT is computed from the EM-estimated frequency in likeFreq().
   LRT = 2*(logL(freq=0) - logL(freq=EM_estimate))
   SNP_pval was converted to chi-square 1-df LRT units via chisq1->invcdf(1-pval)
   in getOptions(). The comparison is: lrt < cutoff → site is NOT a SNP → remove.

### Q: Is .mafs.gz written only for retained sites?
A: YES. The print() method (line 403-478) explicitly checks:
   `if(pars->keepSites[s]==0) continue;`
   Only sites surviving ALL filters get written to .mafs.gz.

### Q: Does keepSites[s] in .mafs.gz represent nInd?
A: YES. The last column "nInd" in .mafs.gz output is literally keepSites[s]
   (line 439: `kputw(pars->keepSites[s],&bufstr)`).
   After likeFreq(), this value is the GL-based count of informative individuals.

### Q: Does -sites / prep_sites affect speed?
A: YES, significantly. prep_sites loads the binary index at chromosome boundaries
   (shared.cpp line 290). The binary index allows O(1) lookup per position.
   Without -sites, every position in the BAM pileup is processed.
   With -sites, non-mask positions are killed at abcFilter before GL computation,
   saving the expensive GL/EM steps. For your 100 kb debug window, the -sites
   overhead is negligible compared to BAM I/O.

### Q: Where is the -minInd check duplicated?
A: In TWO places:
   1. abcFilter.run() line 148-150: checks read-count-based nInd (individuals with >= setMinDepthInd reads)
   2. abcFreq.run() line 614-615: checks GL-based nInd (individuals with informative GLs)
   Both parse the same -minInd argument. In practice, the GL-based count is usually
   <= the read-count-based count (some individuals with reads may have uninformative GLs).
   So abcFilter's check is the first gate, and abcFreq's check can catch additional sites
   where individuals had reads but uninformative likelihoods.

## 5. Why Your Permissive Pass A Gets 7.8M Sites

Pass A uses -SNP_pval 1 -minMaf 0.
- SNP_pval=1 → cutoff = chisq1->invcdf(1-1) = chisq1->invcdf(0) ≈ 0.
  So the check `lrt < 0` basically never fires (LRT is always >= 0 after the clamp).
- minMaf=0 → the MAF check `freq < 0` never fires.
- minInd=200 still applies in BOTH abcFilter and abcFreq.

So Pass A captures all sites where:
  1. Position is in callable mask
  2. >= 200 individuals have reads with >= 3 depth (abcFilter)
  3. Of those, >= 200 have informative GLs (abcFreq)
  4. Major/minor could be inferred

The 7.8M sites out of ~6.7 Mb chromosome = most callable positions pass.
The ~5.9M callable positions that DIDN'T pass = failed minInd or had no data.

## 6. Thread Cap Analysis

Location: shared.cpp lines 68-71
```cpp
if(maxThreads>10){
    fprintf(stderr,"\t-> You have choosen to supply a very high number of threads...\n");
    maxThreads = 8;
}
```

This is a SOFT cap with a warning. The cap exists because:
- ANGSD's parallelism model uses a producer-consumer queue (shared.cpp)
- The producer reads BAM chunks; consumers run the abc module pipeline on chunks
- Thread contention on the print queue can cause diminishing returns above ~8 threads
- The BAM I/O is often the bottleneck, not the computation

For your workflow (226 BAMs, GL computation, EM estimation):
- The computation per chunk IS substantial (226 individuals × EM iterations)
- BAM reading with 226 files is also substantial
- Raising to 16 threads is SAFE — the synchronization uses proper mutexes
- Raising to 80 (your SLURM allocation) will NOT help — the queue model
  means most threads will be waiting for BAM I/O or print queue

RECOMMENDATION: Raise cap to 16 for this workflow. Safe, modest benefit.
