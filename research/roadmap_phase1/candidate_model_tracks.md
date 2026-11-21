# Candidate research gaps and model tracks (Phase 1, Tasks 3, 4 and 6)

**Status:** candidates for the team to evaluate. These are not conclusions. Nothing here has been implemented. Sources are in `references.md`, paper details in `paper_methods_matrix.csv` and datasets in `dataset_feasibility_screen.csv`. ⚠ marks claims that are uncertain or rest on incomplete reading.

---

## Part 1. Research gaps (Task 3)

Five candidate gaps were assessed. **Novelty was checked only against the ~14 microplastic–liver papers reviewed here plus the earlier 19-paper review.** That is not a systematic search, so "nobody has done X" means "not found in these papers".

### G1. Is the hepatic response to plastic diet-dependent? (our original RQ2)
- **Already demonstrated:**
  - Several 2 × 2 studies claim "exacerbation" under HFD/MASH diet: Jung, Fang, Xu P, Ahn, Okamura and Zhai.
  - **Correction to our earlier review:** an interaction *was* computed in Jung (Table S7: 76 genes, mostly skeletal-muscle carry-over), but never reported. Fang and Xu P fitted two-way ANOVAs but printed no interaction statistic. No paper reports and interprets a diet × plastic interaction for molecular data.
  - Our Phase 2 found no robust pathway-level interaction. Gene-level hits were largely muscle artefact.
- **Still uncertain:** whether the absence of an interaction is real or a power failure (median interaction CI ≈ 0.8 log2).
- **Why it matters:** "plastic worsens MASLD" is the field's headline public-health claim, and it rests on a difference-of-significance argument.
- **Relevant methods:** factorial GLMs (DESeq2/edgeR/limma), ROAST/fry and camera with interaction contrasts, equivalence-style reporting (which effect sizes can be excluded).
- **Can our data resolve it?** Only partly. The data can bound the interaction at the program level. They cannot show it is absent.
- **Meaningful new evidence:** program-level interaction estimates with CIs and an explicit "largest interaction compatible with the data", robust to muscle and to leave-one-out (LOO).
- **What makes it inconclusive:** wide CIs at n = 4–5, with no second 2 × 2 microplastic dataset with public omics anywhere (Fang, Xu P and Ahn deposited nothing).
- **Verdict:** mostly done already. The remaining increment is in reporting (bounds, equivalence), not a new model.

### G2. Does the PE-associated signature resemble MASH, or another physiological state such as PPARα activation or fasting/feeding/circadian phase? (our original RQ1, sharpened)
- **Already demonstrated:**
  - Jung shows overlap between PE and MD DEG lists, using pairwise comparisons against CD+Veh.
  - Our split-reference analysis found a modest PE(CD)–diet resemblance (ρ ≈ 0.33) that is absent under MD.
  - That resemblance is built from the **same mice**, and the shared CD+Veh reference has a distinct clock state.
  - The PE program pattern is lipogenesis and SREBP-cholesterol down with PPARα targets up. That is also the textbook fasting signature.
- **Still uncertain:** which reference state the PE signature is actually closest to, once the disease axis comes from **independent** mice.
- **Why it matters:**
  - "PE mimics MASH" and "PE induces a fasting-like or PPARα state" imply different biology: disease-like injury versus altered feeding or energy status, possibly gut-mediated.
  - Several studies suggest an indirect (gut) route: Zhai (FMT), Park (no PET found in liver), Djouina.
- **Relevant methods:**
  - signature projection and direction-agreement scoring against reference contrasts (Vacca 2024, Teufel 2016);
  - multiple regression of one contrast's fold changes on several reference contrasts, with model selection.
  - ⚠ Connectivity-map-style projection is standard in toxicogenomics, but its primary sources were not reviewed here.
- **Can our data resolve it?** **Partly, and better than before.**
  - GSE316023 uses the **same MASH diet (D09100310)**, giving an independent MASH axis (GAN vs chow, 8 vs 8 at week 23).
  - GSE246328 is a second, independent GAN axis that includes an 8-week time point.
  - GSE282998/GSE8295 give a PPARα activation axis.
  - A fasting/feeding or time-of-day reference **was not screened yet** ⚠. It must be found in Phase 2, or this sub-question drops.
- **Meaningful new evidence:**
  - On independent axes, the PE contrast is better explained by one reference axis than by others.
  - That explanation holds under both diets, and passes a positive control: our own Diet_in_Veh contrast must project strongly onto the external MASH axis.
- **What makes it infeasible or inconclusive:**
  - The PE signature itself is confounded with library depth and yield. Projection inherits this; it can only say what the signature looks like, not that PE caused it.
  - No fasting/circadian reference may exist with counts.
  - Chow (not purified-diet) comparators in the MASH references inflate the diet axis.
  - Overlapping axes (PPARα and fasting share targets) may make the model choice unidentifiable.
- **Verdict:** **strongest candidate.** It directly upgrades RQ1, uses independent data, and maps onto course regression and model selection.

### G3. Are microplastic liver responses reproducible across laboratories at the program level, and is the lipogenesis direction conflict systematic?
- **Already demonstrated:** lipogenesis directions conflict.
  - **Down:** Jung (PE), Makdissi (PS, module level), Chiu (PS, abstract).
  - **Up:** Djouina (PE, qPCR), Okamura (PS + HFD, qPCR), Xu P (PS-NP), Li (Fas, by qPCR).
  - PPARα is up in Li and Djouina, suppressed in Fang (after an overnight fast), and not significant within MD in Jung.
  - Our 8-gene replication found Djouina's lipogenic increases opposite in direction.
- **Still uncertain:** whether the conflict is driven by assay (qPCR panel vs RNA-seq), sampling or fasting state, diet, route, polymer or size. Or whether there is simply no reproducible signal.
- **Why it matters:** without cross-lab reproducibility, no "microplastic liver signature" can be claimed (earlier gap C5).
- **Relevant methods:** within-study effect estimation followed by synthesis: sign and vote counting, rank-based methods (RankProd), Hedges' g random effects (Ramasamy 2008; Haynes/Sweeney 2017), leave-one-study-out checks.
- **Can our data resolve it?** **Only descriptively.** Independent public liver count data with a plastic contrast:
  - GSE245069 (PS, ICR, 4 vs 4 × 3 sizes, one shared control);
  - GSE277772 (PS, 3 vs 3; counts availability conflicting ⚠);
  - the Hong-lab series GSE261893/GSE261725/GSE286116/GSE262399 (PE/PS/PP/PB, 3 each, **byte-identical shared controls**, 71 weeks, no paper).

  That is roughly **4 independent laboratories**, each with 3–5 mice per arm. Polymer, strain, diet, route and duration all vary together across labs.
- **Meaningful new evidence:**
  - Pre-specified programs (lipogenesis/SREBP1, SREBP2-cholesterol, PPARα targets, inflammation) move in a consistent direction across labs.
  - Or the evidence shows they do **not**, which is equally publishable.
  - Plus a demonstration that program-level agreement exceeds gene-level DEG overlap (as Teufel and Vacca found for mouse vs human).
- **What makes it inconclusive:**
  - n = 3 per arm in two of the four labs;
  - between-study variance cannot be estimated with ~4 studies (Sweeney 2017 never tested this regime);
  - moderators are confounded, so heterogeneity cannot be attributed to polymer or diet.
- **Verdict:** feasible as a **descriptive reproducibility study**. Not feasible as pooled inference.

### G4. Are technical library differences aligned with treatment in other microplastic transcriptomes?
- **Already demonstrated:** in GSE294180, PE completely separates from vehicle on depth and yield (our audit). No reviewed paper reports sample-level RNA-seq QC (earlier review, gap A3).
- **Still uncertain:** whether this is specific to Jung's dataset or common in the field.
- **Why it matters:** if common, published "plastic signatures" may partly be library artefacts.
- **Methods:** read-depth and yield metadata from SRA/ENA, regressed or tested against group; library-composition checks.
- **Feasible?** Yes, cheaply: run metadata exist for every dataset in the screen.
- **What counts as evidence:** treatment–depth association statistics per study.
- **Inconclusive if:** an association exists but its effect on DE cannot be quantified (the same identifiability problem we already face).
- **Verdict:** a small, solid **add-on QC step for Track B**, not a project on its own.

### G5. Composition versus per-cell expression
- **Demonstrated:** Makdissi shows Kupffer-cell mechanisms, and Jung shows spatial shifts in cell proportions.
- **Our data:** bulk only. The single-nucleus reference GSE291616 has n = 2, and the Xenium deposit lacks mouse IDs and may confound slide with PE (unconfirmed).
- **Verdict:** **out of scope for Phase 2** unless the team wants a deconvolution side analysis. Record it as a limitation.

---

## Part 2. Candidate model tracks (Task 4)

For each track: **Published paper → existing method → demonstrated limitation → proposed improvement → appropriate data → baseline → evaluation metric → external benchmark.**

### Track A: Factorial, program-level decomposition of the PE signature against independent reference axes (G2 + G1)

**Question.** Within GSE294180, is the PE-associated hepatic response better described as MASH-like, PPARα-activation-like or fasting/circadian-like? And does that description differ between diets?

| Link | Content |
|---|---|
| Published paper | Jung 2026 (PE "compounds" MASH; DEG overlap with the MD signature) |
| Existing method | Pairwise DESeq2 vs CD+Veh, DEG set subtraction, WGCNA; unreported interaction table (S7) |
| Demonstrated limitation | Pairwise contrasts mix diet and PE. The disease axis comes from the same mice, so shared controls bias the PE–MD overlap (our split-reference analysis showed ≈ 40% of the naive correlation is artefact). The Ppara headline is vs CD+Veh. S7 is muscle-dominated. |
| Proposed improvement | (1) Keep the factorial model and gene-set tests (fry primary, camera supporting) for program-level PE effects under each diet, and report interaction CIs as bounds. (2) Build **external reference axes**: MASH (GSE316023 GAN vs chow wk23; GSE246328 GAN 8 wk vs chow), PPARα activation (GSE282998), fasting/circadian (to be found ⚠). (3) Regress PE_in_CD and PE_in_MD fold changes on the reference fold changes (weighted by SE). Compare single-axis and multi-axis models by AIC and by **gene-set-blocked cross-validation**, so that held-out pathways are predicted. (4) Stratified bootstrap over mice for CIs. |
| Appropriate data | GSE294180 (deposited counts, assumptions A1–A3) plus the reference datasets above |
| Baseline | (a) Jung-style DEG overlap with MD; (b) our same-data split-reference ρ = 0.33 |
| Evaluation metric | Variance in PE fold change explained by each axis (held-out gene sets); coefficient CIs; for G1, interaction estimate and CI per program |
| External benchmark | **Positive control:** our Diet_in_Veh contrast must align strongly with the external MASH axes, otherwise the axis is unusable. **Negative and technical controls:** binomial thinning of vehicle libraries to PE depth, to test whether depth alone reproduces the projection; and a muscle-score–adjusted fit. |

- **Strengths:** removes the same-data circularity. Turns RQ1 into a competing-explanations question that addresses our biggest alternative explanation (clock/feeding). Directly uses regression, model selection, CV, bootstrap and hypothesis testing. Cheap: counts are public.
- **Weaknesses:** cannot break the PE ≡ depth confound. Thinning tests depth but not library *quality*. Reference axes come from different labs, strains and comparators (chow). PPARα and fasting axes overlap.
- **Kill criteria:** the external MASH axis fails the positive control; or no fasting/circadian reference with counts can be found (then the question narrows to MASH vs PPARα); or axes are too collinear to compare.

### Track B: Cross-study program-level reproducibility of microplastic liver responses (G3 + G4)

**Question.** Do lipogenic, cholesterol-synthesis, PPARα and inflammatory programs move consistently with plastic exposure across independent laboratories, and is program-level agreement better than gene-level DEG overlap?

| Link | Content |
|---|---|
| Published paper | Li 2025 (CytoHubba hub genes from DEG overlap); Djouina 2023 (qPCR panel); our 8-gene replication |
| Existing method | DEG lists at study-specific thresholds; Venn overlap; PPI hub ranking; single-gene qPCR comparison |
| Demonstrated limitation | Thresholds differ by orders of magnitude between papers. Hub status is topology, not evidence (Fmo2 failed qPCR). Gene-level directions conflict across labs. DEG overlap is known to under-state agreement (Teufel 2016: 1–18 of 65–177 genes shared; pathways agree better). |
| Proposed improvement | Uniform re-processing from counts (or SRA where counts are missing); within-study DESeq2 with each study's own controls; pre-specified program scores (fry/camera effect plus mean log2FC); synthesis by sign concordance, rank-based statistics and a descriptive random-effects estimate with heterogeneity reported, **not** used for inference; leave-one-lab-out; SRA depth-vs-group QC per study (G4). The Hong-lab series are treated as **one study** with shared controls, and GSE245069's three sizes as one study with correlated contrasts. |
| Appropriate data | GSE294180 (PE_in_CD, PE_in_MD), GSE245069, GSE277772 (verify counts), Hong-lab GSE261893 + GSE261725 + GSE286116 + GSE262399. **Do not pool raw counts.** ComBat-seq is invalid when lab is nested with treatment (Zhang 2020; Leek 2010). |
| Baseline | Gene-level DEG overlap (Jaccard) and hub-gene agreement between studies |
| Evaluation metric | Pairwise between-study agreement (sign agreement on pre-specified programs; Spearman ρ of program effects) for program level vs gene level; number of programs with consistent direction across ≥ 3 labs |
| External benchmark | No gold standard exists. The nearest is a positive-control replication: diet-axis agreement between GSE316023 and GSE246328 (two GAN studies) gives an upper reference for what "reproducible" looks like with similar methods |

- **Strengths:** speaks to the field's actual problem (conflicting directions). Uses only public data. A negative result is still informative.
- **Weaknesses:** about 4 labs, with n = 3–5. Heterogeneity is uninterpretable because polymer, strain, diet, route and duration all co-vary. PE is represented by only 2 labs, one with no paper and unknown particle size.
- **Kill criteria:** GSE277772 counts cannot be obtained and SRA re-quantification is out of scope (then 3 labs); program scores in n = 3 studies have CIs spanning everything.

### Track C: Predictive machine learning, assessed and mostly rejected

**Options considered and rejected**

| Target | Why rejected |
|---|---|
| PE vs vehicle classifier within GSE294180 | Perfect separation by sequencing depth and yield; any classifier can learn library properties. With n = 19 there is no held-out set, and CV with feature selection is optimistically biased (Ambroise & McLachlan 2002; Varma & Simon 2006). |
| "Microplastic-exposed" classifier pooled across studies | Each lab contributes one treatment pattern, so exposure is nested in lab. A classifier would learn the lab, polymer or platform. Cross-study validation overstates transfer even with 100–500 samples per study (Bernau 2014); here studies have 6–19 samples. |
| Polymer or particle-size prediction | One lab per polymer–size combination; not identifiable |
| Toxicogenomic hepatotoxicity classifier | Published successful examples needed hundreds of arrays (Fielden 2007; Natsoulis 2005). There are no labelled hepatotoxicity outcomes per mouse in our data. |

**The only scientifically meaningful option: an externally trained MASH-state score applied to GSE294180**

| Link | Content |
|---|---|
| Published paper | Jung 2026 (PE "MASH-like"); Vacca 2024 (ranking mouse models by proximity to human MASLD) |
| Existing method | DEG overlap; direction agreement with hypergeometric tests |
| Demonstrated limitation | Same-data disease axis; no out-of-sample test |
| Proposed improvement | Train a regularised logistic regression (elastic net) or a simple signed-gene score on GAN vs chow in GSE316023, with nested CV. Validate on GSE246328 (cross-study). Then score GSE294180: Diet_in_Veh mice are the **positive control** (they should score MASH-like), and PE vs vehicle within each diet is the question. |
| Appropriate data | Training: GSE316023 (GAN 62, chow 16; chow only at week 23, so restrict to wk23 or model time explicitly). Validation: GSE246328. Test: GSE294180. |
| Baseline | Unweighted mean z-score of the top external diet genes (no ML) |
| Evaluation metric | Cross-study AUC and calibration on GSE246328; AUC for MD vs CD vehicle in GSE294180; score difference PE − Veh with a bootstrap CI |
| External benchmark | GSE294180 diet labels (an independent lab, the same diet as GSE316023) |

- **Assessment:** this is **Track A's MASH axis implemented as a classifier**, not an independent direction. Unless elastic net beats the simple score on cross-study validation, the simpler score should be used.
- **Recommendation:** fold it into Track A as an alternative implementation and sensitivity analysis. **Do not run it as a separate track.**
- **Remaining risks:** the chow comparator inflates the axis; class imbalance and a time confound in GSE316023; scored PE samples are still depth-confounded.

---

## Part 3. Where the BT3041 course methods genuinely contribute (Task 6)

| Course technique | Where it earns its place | Where it should *not* be forced |
|---|---|---|
| Hypothesis testing and multiple testing | Factorial contrasts; program-level tests (fry/camera); positive and negative controls (A); sign tests across labs (B) | Gene-by-gene cross-study "significance" counting |
| Regression | NB-GLM factorial model; **regression of PE fold changes on reference axes** (A, core); additive vs factorial comparison | Regression with depth as a "correction" (it is collinear; use only as a stress test) |
| Model selection | AIC/LRT for single vs multi-axis explanations (A); additive vs interaction | Picking among dozens of ML models on n = 19 |
| Bootstrapping | Stratified mouse bootstrap for axis coefficients and program effects (A); leave-one-lab-out (B) | Treating the bootstrap as a replacement for more mice (percentile CIs are too narrow at n = 4–5) |
| Cross-validation | Gene-set-blocked CV for axis models (A); nested CV for the MASH score trained on external data (C, inside A) | CV within GSE294180 for exposure classification (rejected) |
| Classification | Only the externally trained MASH-state score, with cross-study validation | PE-exposure or polymer classifiers |
| Clustering | QC and contamination detection (already done); grouping reference axes or programs descriptively | Discovering "PE subtypes" in 19 mice |
| PCA | QC (done); summarising program scores across studies (B) | Drawing evidence from PCA of the same mice used for testing |

---

## Part 4. Summary comparison

| | Track A | Track B | Track C |
|---|---|---|---|
| Scientific gap | G2 (+ G1 reporting) | G3 (+ G4 QC) | Part of G2 |
| Data in hand or public | Yes (fasting reference still to find) | Yes (GSE277772 counts to verify) | Yes |
| Independent evidence | Yes (external axes) | Yes (other labs) | Yes (external training) |
| Breaks PE ≡ depth confound? | No (thinning addresses depth only) | Partly: other labs have their own artefacts, so a direction consistent across labs is harder to explain by one lab's artefact | No |
| Main risk | Axis collinearity; missing fasting reference | Too few, too small, too heterogeneous studies | Duplicates A; overfitting |
| Course fit | Strong (regression, model selection, bootstrap, CV) | Moderate (testing, bootstrap, PCA) | Classification/CV, but the least scientifically necessary |
| Suggested status | **Primary candidate** | **Secondary or complementary candidate** | **Merge into A; do not run separately** |
