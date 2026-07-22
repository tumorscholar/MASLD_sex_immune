# The statistics behind the sex × MASLD meta-analysis — explained plainly

*A walk-through of every method in the pipeline: what it does, why it is the right tool for this kind of biological data, and where its limits are. Written to be read start-to-finish, no statistics background assumed.*

---

## The one-paragraph version

We took five public liver datasets, worked out each patient's sex from their gene expression, turned each immune gene-signature into a single score per patient, and then asked one question for every signature: **once we account for fibrosis stage and for the fact that the datasets are different, is this immune readout systematically higher in men or in women?** Everything else in the pipeline exists to make that question answerable honestly — to stop different datasets from contaminating each other, to stop one study driving the whole result, to stop metabolic differences masquerading as sex, and to check the answer with a second, independent method and a disease-free comparison. None of the methods are exotic. They are the standard, defensible toolkit for pooling transcriptomic studies, and I have tried to explain below not just *what* each one is but *why* it is the appropriate choice for data like ours.

---

## Step 0 — What the data actually are, and why that dictates the statistics

Each dataset is a big table: genes down the side, patients across the top, and in each cell a number telling us how strongly that gene was switched on in that patient's liver. Two of our cohorts are RNA-seq (the number is a read count) and one is a microarray (the number is a fluorescence intensity). This matters for the whole analysis, because it forces two rules that everything downstream obeys:

1. **The absolute numbers are not comparable between datasets.** A count of 500 in an RNA-seq study and an intensity of 500 on a microarray mean completely different things. So we can never pool raw values — we can only pool *relative* patterns within each dataset. This is why you will see z-scoring and rank-based scoring everywhere: they both throw away the arbitrary scale and keep the biology.

2. **Patients from the same study are more alike than patients from different studies**, for reasons that have nothing to do with sex or fibrosis — the same hospital, the same sequencing machine, the same batch of reagents. Statistically these patients are *not* independent, and pretending they are would make our p-values far too optimistic. This is the single most important reason we use a mixed-effects model rather than a plain regression (explained in Step 4).

Everything below follows from these two facts.

---

## Step 1 — Working out each patient's sex from their genes (`01_sex_assignment`)

**What it does.** Rather than trusting the sex label in the metadata (which is sometimes missing, sometimes wrong), we recover sex directly from the expression data. We look at *XIST* — a gene switched on almost exclusively in cells with two X chromosomes — and a panel of nine Y-chromosome genes (*RPS4Y1, DDX3Y, EIF1AY, UTY, KDM5D, USP9Y, NLGN4Y, ZFY, TXLNGY*). Women have high XIST and silent Y-genes; men are the mirror image.

**The statistic.** Within each cohort we *z-score* XIST and the averaged Y-gene signal. A z-score just re-expresses each value as "how many standard deviations above or below this cohort's average" — so +1 means one SD above the mean, −1 one SD below. We then take the difference (Y minus XIST). A clearly positive difference is male, a clearly negative one female, and a small band in the middle (within 0.5 SD) is left "ambiguous" and dropped rather than guessed.

**Why this is the right tool for biological data.** Sex is one of the very few things we can verify against an external truth in transcriptomic data, because the X/Y genes are biological ground truth, not a soft signature. Z-scoring *within* each cohort is what makes it work across a microarray and two RNA-seq platforms at once — because a z-score is unit-free, a male on the array and a male on RNA-seq both land in the same place. Leaving the ambiguous middle band unassigned is deliberate conservatism: a few uncertain calls do no harm, but a confidently *wrong* sex call would.

**What to keep in mind.** This is an expression-based *proxy* for chromosomal sex; it can't see sex-chromosome aneuploidies (e.g. XXY) and it very occasionally mislabels a low-quality sample. Where the metadata *did* record sex, we checked our calls against it and agreement was high — that concordance is our validation that the method is trustworthy.

---

## Step 2 — Turning a gene list into one number per patient: singscore (`02_build_matrix`)

**The problem.** "Are MAIT cells more common in men?" is a question about a *cell type*, but bulk liver data doesn't count cells — it measures genes. So for each cell type or functional state we define a small marker-gene set (e.g. MAIT ≈ *SLC4A10, KLRB1, ZBTB16, RORC, TRAV1-2*) and need to collapse those genes into a single "how much MAIT signal does this patient have" number.

**The statistic — singscore, in plain terms.** Within one patient, rank every gene from lowest-expressed to highest. A marker gene that sits near the top of that patient's ranking is strongly on; near the bottom, off. The score is essentially **the average rank of the marker genes, rescaled to run from −0.5 (all markers at the bottom) through 0 (markers sitting exactly mid-pack) to +0.5 (all markers at the top).** For signatures that have both "up" and "down" genes (like tissue-resident memory T cells, which turn some genes on and others off), we score the up-set and the down-set separately and subtract.

**Why rank-based scoring is the appropriate choice.** Because it uses *ranks*, not raw values, it is almost immune to the two problems from Step 0 — different platforms, different normalisations, different dynamic ranges all wash out when you only care about "is this gene high or low *for this patient*." It needs no training data and no reference cohort; each patient is scored entirely on their own profile, so there is no way for one dataset to leak into another's scores. That self-contained, non-parametric quality is exactly what you want when you are going to pool heterogeneous studies. It is a published, widely-used method (Foroutan et al., 2018) precisely for single-sample scoring of this kind.

**What to keep in mind.** A signature score is a *relative* quantity — it tells you a patient is high or low compared to the rest of that cohort, not an absolute cell fraction. That is fine for our question (which is comparative) but it is why we never report "X% MAIT cells" from the bulk data. It also means a score is only as good as its marker set; we chose canonical, literature-standard markers and later cross-checked the key results with a completely different method (Step 7).

---

## Step 3 — Making cohorts comparable before pooling (within-cohort z-scoring of scores)

Before any modelling, every signature score is z-scored *again, within its own cohort*. This is the step that lets us pool five studies without one study's baseline swamping the others. After this, a score of "+1" means the same thing — one SD above that study's average — in every cohort. We are deliberately comparing *relative position within a study*, never absolute values across studies. This is the meta-analytic equivalent of grading each class on its own curve before comparing students across schools.

---

## Step 4 — The core test: the linear mixed-effects model (`03_maineffect`)

This is the heart of the paper, so it's worth going slowly.

**The question in one line:** holding fibrosis stage constant, is each immune readout different between men and women, across all cohorts together?

**The model:** `readout ~ sex + fibrosis + (1 | cohort)`

Read it left to right. We are explaining the immune readout using two things we care about — **sex** and **fibrosis stage** — plus a bookkeeping term for **which cohort** the patient came from.

- **`sex` and `fibrosis` are "fixed effects."** These are the specific, repeatable things we want to measure. The number the model gives back for sex (called *beta*, β) is the answer: how much the readout shifts, in standard-deviation units, going from female to male, *after fibrosis has been accounted for*. A positive β means higher in men, negative means higher in women. Putting fibrosis in the same model is what makes the sex estimate "stage-adjusted" — we are not confusing a sex effect with the fact that one sex might be sicker.

- **`(1 | cohort)` is a "random effect" — the crucial part.** It tells the model: "each cohort is allowed its own baseline level; don't treat the five studies as one big pool of independent patients." This directly fixes problem 2 from Step 0. Without it, two patients from the same study would be counted as two independent pieces of evidence when they are really partly redundant, and the p-values would be dishonestly small. The random intercept absorbs that shared, study-level baseline so that the sex estimate is driven by *within-cohort* male-vs-female differences, then combined across cohorts.

**Why a mixed model is the correct tool here, not a t-test or plain regression.** A t-test would ignore fibrosis entirely and pretend all 481 samples are independent — both wrong. A plain (fixed-effect) regression could adjust for fibrosis but still treats cohorts as interchangeable. The mixed-effects model is the standard, textbook way to analyse data that come in *groups* (here, studies) where you want a general effect while respecting that the groups differ. It is exactly the structure used in formal gene-expression meta-analyses. In effect it does an internal meta-analysis: estimate the sex effect within each cohort, then pool, weighting appropriately.

**What the model hands back for each readout:** the effect size **β** (direction and magnitude, in SD units), its **standard error** (how precise that estimate is), and a **p-value** (how surprised we'd be to see an effect this large if sex truly made no difference).

**Two fibrosis codings, on purpose.** We run the whole thing twice — once treating fibrosis as an ordered scale (F0→F4) and once as a simple yes/no (any fibrosis vs none). If a result only appears under one arbitrary coding choice, it's fragile. We only call a readout a "headline" finding if it survives *both* codings with the *same direction*. That is a built-in sensitivity analysis, not extra decoration.

---

## Step 5 — Guarding against false positives across many tests: BH-FDR

We don't test one readout, we test ~15. If you run 15 independent tests at the usual p < 0.05, you expect roughly one "significant" hit by pure chance even if nothing is real. Reporting that as a discovery would be a mistake.

**The fix — Benjamini–Hochberg false discovery rate (FDR).** Instead of controlling the chance of *any* false positive (which is very strict and costs you real findings), BH controls the *proportion* of your declared hits that are false. An FDR of 0.10 means "of the readouts I'm calling significant, I expect about 10% to be flukes." For discovery-stage biology — where we want a shortlist of promising, sex-biased readouts to then validate — FDR is the appropriate and widely-accepted middle ground between missing everything (Bonferroni) and believing everything (raw p-values). Every headline number in the paper is an FDR-adjusted value, not a raw p.

---

## Step 6 — Is one study secretly driving the result? Leave-one-cohort-out (LOCO)

A pooled result can be an illusion created by a single unusual dataset. So for every readout we refit the model five times, each time **dropping one whole cohort**, and check two things: does the direction of the sex effect stay the same, and does it stay significant. A finding that flips direction or evaporates when you remove one study is not robust; a finding that holds no matter which study you leave out is. This is a leave-one-out robustness check, and it is one of the most honest things you can do to a meta-analysis — it directly answers the reviewer's inevitable "is this just cohort X?"

---

## Step 7 — A second opinion from a different method: deconvolution concordance (`06_prep_deconv`, `07_run_deconv`, `08_concordance`)

Singscore is one way to read cell content out of bulk tissue. To make sure our sex directions aren't an artefact of that one method, we re-estimate cell content with two **completely independent, published deconvolution tools** — **xCell** and **MCP-counter** — which use their own reference signatures and their own maths, then run the *same* mixed model on their outputs.

**The logic — concordance, not agreement in magnitude.** We are not asking the three methods to produce the same number; different methods never do. We are asking whether they **agree on the direction** of each sex effect — if singscore says "Tregs higher in women," do xCell and MCP-counter independently point the same way? Agreement across methods that share no code and no assumptions is strong evidence the biology is real and not a scoring quirk. We report the fraction of matched cell types where the directions agree.

**Why this is a fair test.** xCell and MCP-counter were trained largely on blood and tumour, not liver, so they are genuinely independent of our marker choices — they can't just be echoing our singscore signatures. That independence is the whole point. It's also why one cell type, MAIT, *can't* be checked this way — it's defined by its T-cell receptor and simply isn't in these deconvolution references, so we validate MAIT differently (receptor-identity genes, the GTEx result, and single-cell data) and say so plainly.

---

## Step 8 — Is it really sex, or is it BMI and diabetes? Deconfounding (`09_deconfound`)

Men and women in these cohorts differ in more than sex — often in BMI and type-2-diabetes rates too. If we're not careful, a "sex effect" could actually be a metabolic effect wearing a sex costume.

**The test.** In the one cohort that records BMI, T2D and age (GSE89632), we fit the sex effect twice on the *exact same patients*: once unadjusted, once adding BMI, T2D and age to the model. Then we look at how much the sex β shrinks. If it barely moves, metabolic state is *not* confounding the sex effect. Doing this within a single cohort is deliberate — it isolates the effect of *adjustment* from the effect of *losing samples*, so a change can only mean confounding, not lost power.

**Why this is the right design.** The clean way to ask "does B explain away A?" is to hold the sample fixed and add B to the model — anything else confuses confounding with a change in who's being analysed. This is a standard deconfounding / covariate-adjustment argument.

---

## Step 9 — Is the sex difference constitutional or disease-driven? The GTEx control (`04_gtex_control`)

This is the conceptual keystone of the paper's two-axis story. For a sex difference in MASLD liver, there are two very different explanations:

- **Constitutional:** men and women differ in this immune readout in *healthy* liver too, and MASLD just inherits it.
- **Disease-emergent:** the sex difference only appears once disease is present.

**The test.** We score the *identical* immune readouts in **GTEx** — hundreds of disease-free liver samples — and ask, for each headline readout, whether the healthy liver shows the *same* sex direction (→ constitutional) or *no/opposite* difference (→ MASLD-specific). Because we use the very same signatures and the same kind of model (sex + age), it's an apples-to-apples comparison against a disease-free baseline.

**Why this matters and where it's fragile.** This is what lets us say a male-MAIT bias is a baseline trait the liver carries into disease, whereas a female Treg/CD8/DC signal *emerges with* disease. It is a genuinely strong design. Its honest weakness — and the one the QED reviewer rightly pressed — is that GTEx liver is *post-mortem*, immune-cell-poor, and not sex-balanced or adjusted for tissue-quality variables. That's exactly why the agreed next step is to re-run the GTEx comparison with technical covariates (RIN, ischaemic time, centre) and sex-balancing, so the disease-free baseline is as clean as the disease analysis.

---

## Step 10 — How to read the numbers in the paper

- **β (beta), "higher in males/females":** the effect size in standard-deviation units, after adjustment. Sign = direction; size = how big. β = 0.4 means the readout sits about 0.4 SD higher in one sex — a modest but real shift, typical for immune signatures in bulk tissue.
- **p-value:** probability of seeing an effect this large by chance if sex truly did nothing. Smaller = more surprising-if-null.
- **FDR:** the p-value *after* correcting for testing many readouts; this is the one to trust for the headline claims.
- **LOCO-consistent:** the direction held when each cohort was left out in turn — i.e. not driven by one study.
- **Concordant (deconvolution):** an independent method agreed on the direction.

A finding earns "headline" status only when it clears **all** of these at once — both fibrosis codings, FDR, LOCO, and (where possible) a second method. That is a deliberately high bar, and it's why the paper's central claims are few but defensible.

---

## The honest summary on "are these methods appropriate for biological data?"

Yes — and deliberately so. Every choice in the pipeline is the conservative, standard-practice option for its job: rank-based scoring because it survives platform differences; z-scoring because it makes cohorts comparable; mixed-effects models because the data come in study-groups and independence is false; FDR because we test many readouts; LOCO and dual codings because single studies and arbitrary codings shouldn't decide a result; independent deconvolution because one method shouldn't either; and a disease-free control because "different in disease" only means something against a healthy baseline. The methods are not the weak point of this paper. The weak points are the ones we name ourselves — bulk deconvolution isn't liver-trained, signature scores are relative not absolute, expression-based sex is a proxy, and the GTEx baseline needs its technical covariates — and each of those has a planned, concrete fix. That is the right way round: strong, boring, defensible statistics, with the caveats stated out loud rather than hidden.
