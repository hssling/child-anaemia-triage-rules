# Protocol v1.0 — Household triage rules for child haemoglobin testing in India

**Manuscript ID:** 12_child_anaemia_triage_rules
**Programme:** Anaemia in India Research Programme
**Protocol version:** 1.0
**Date frozen:** 5 August 2026
**Freezing rule:** this file is SHA-256 hashed and the hash recorded in
`governance/PROTOCOL_FREEZE.txt` *before* the confirmatory analysis script is
executed. Any change after that point creates v1.1 and must be declared as a
post-hoc amendment in the manuscript.

---

## 1. Rationale

India's Anemia Mukt Bharat (AMB) framework was revised in June 2026 into a
7×7×7 structure built on T4 — Test, Treat, Talk, Track. The first T is Test.

Two facts create the question this protocol answers.

1. The AMB health management information system contains **no haemoglobin
   testing indicator for children aged 6–59 months**. The only two indicators
   for this beneficiary group are IFA syrup doses and deworming. Haemoglobin
   measurement indicators exist solely for pregnant women. A pathway whose
   first step is Test does not count the test in the largest child beneficiary
   group.
2. The target population is large — approximately 111.8 million children aged
   6–59 months — and point-of-care haemoglobin testing consumes a consumable
   cuvette, a lancet and front-line worker time per child.

Where universal testing is not delivered, front-line practice necessarily
falls back on informal prioritisation. The mother's anaemia status is the most
commonly available such signal: mothers are tested under the maternal arm of
the same programme, and maternal and child anaemia are known to be
statistically associated.

Association, however, is not a decision rule. No published study has converted
the maternal–child anaemia association into the quantities that a triage
decision actually requires: sensitivity, negative predictive value, likelihood
ratios and net benefit. This protocol specifies that evaluation.

## 2. Objective

**Primary objective.** To determine whether maternal anaemia status, alone or
in simple pre-specified combinations with routinely available household
characteristics, is an acceptable triage rule for deciding which children aged
6–59 months receive a haemoglobin test in India.

**Secondary objectives.**
1. To quantify, for each candidate rule, the trade-off between tests avoided
   and anaemic children missed.
2. To determine whether any Indian state or union territory has a maternal-
   status rule that meets the acceptability bar even if India as a whole does
   not.
3. To estimate the net benefit of each rule relative to testing all children
   and testing no children across the plausible range of decision thresholds.

## 3. Design

Cross-sectional diagnostic-accuracy and decision-analytic evaluation of
pre-specified triage rules in a nationally representative survey, reported
according to STARD 2015 with an explicit decision-curve extension.

This is a **rule-evaluation** study, not a determinants study. No causal
relationship between maternal and child anaemia is estimated, asserted or
implied. No inference about iron deficiency is made from haemoglobin.

## 4. Data source

National Family Health Survey 2019–21 (NFHS-5), India, Children's Recode
(`IAKR7EFL.DTA`), obtained from The DHS Program under an approved data request.
Held locally at `NFHS 5 microdata/IAKR7EDT.zip`.

NFHS-5 measured haemoglobin using the HemoCue Hb 201+ point-of-care system on
capillary blood. All haemoglobin values used here are the survey's
altitude-adjusted values; maternal values are additionally adjusted for
smoking, per DHS standard recoding.

## 5. Participants

**Inclusion.** Living children aged 6–59 completed months at measurement, with
a valid altitude-adjusted haemoglobin value, linked to their interviewed mother
who also has a valid adjusted haemoglobin value.

**Exclusion.** Children not measured, measured with an invalid or missing
result, outside the 6–59 month window, or whose mother has no valid
haemoglobin value.

**Expected analytic sample.** Approximately 182,000 mother–child pairs.

**Unit of analysis.** The child. Mothers may contribute more than one child;
the variance estimator accounts for clustering at the primary sampling unit,
which subsumes clustering within household and within mother.

## 6. Reference standard (target condition)

Child anaemia, defined as altitude-adjusted haemoglobin **< 11.0 g/dL** in
children aged 6–59 months, the threshold used by WHO and by the AMB programme.

Severity strata for descriptive reporting only: mild 10.0–10.9, moderate
7.0–9.9, severe < 7.0 g/dL.

The reference standard is measured in every participant, so there is no partial
or differential verification.

## 7. Index tests — pre-specified triage rules

Each rule is a binary decision to test the child. All inputs are facts a
front-line worker plausibly holds without performing the child's test.

| Rule | Definition (test the child if…) |
|---|---|
| R1 | Mother is anaemic |
| R2 | Mother is moderately or severely anaemic |
| R3 | Child is aged 6–23 months |
| R4 | Mother is anaemic **OR** child is aged 6–23 months |
| R5 | Mother is anaemic **OR** household is in the poorest two wealth quintiles |
| R6 | Mother is anaemic **OR** child is aged 6–23 months **OR** household is in the poorest two wealth quintiles |
| R7 | Maternal haemoglobin below a threshold, evaluated continuously (ROC) |

Maternal anaemia is defined with the pregnancy-appropriate threshold: adjusted
haemoglobin < 11.0 g/dL if currently pregnant, < 12.0 g/dL otherwise.
Maternal moderate-or-severe anaemia is < 10.0 g/dL if pregnant, < 9.0 g/dL
otherwise.

Rules R1–R6 are fixed in advance. R7 is exploratory and is used only to
establish the best achievable performance of *any* maternal-haemoglobin
threshold, so that the conclusion cannot be attributed to a poor choice of
cut-off.

## 8. Acceptability bar — declared before results are read

A triage rule is declared **usable** only if it achieves, at national level,
both of:

- **Sensitivity ≥ 90%** — no more than one in ten anaemic children is missed; and
- **Negative predictive value ≥ 80%** — of the children a rule declines to
  test, at least four in five are genuinely non-anaemic.

Both criteria must be met by the point estimate **and** the lower bound of the
95% confidence interval must exceed 85% and 75% respectively.

This bar is stated here, before analysis, so that it cannot be relaxed after
seeing the results.

## 9. Statistical methods

**Survey design.** All estimates use the NFHS-5 complex design: individual
women's weight `v005`/1,000,000, primary sampling unit `v021`, and sampling
stratum `v022`. Finite population correction is not applied.

**Point estimates.** All accuracy measures are ratios of weighted totals.
Sensitivity, specificity, PPV and NPV are estimated as design-weighted ratio
estimators.

**Variance.** Design-based variances are obtained by Taylor linearisation of
the ratio estimator. For a ratio R = X/Y of weighted totals, the linearised
variate is z = w(x − R·y), and the variance of its total is computed with the
stratified, ultimate-cluster formula

  Var(T) = Σ_h [ n_h/(n_h − 1) · Σ_i (t_hi − t̄_h)² ]

where t_hi is the PSU-level total of z and n_h is the number of PSUs in stratum
h. Strata containing a single PSU are collapsed with the nearest stratum in the
same state and their treatment is reported.

Confidence intervals for proportions are computed on the logit scale and
back-transformed, to keep intervals inside [0, 1].

**ROC and AUC.** The weighted ROC curve for continuous maternal haemoglobin is
computed by design-weighted cumulative sensitivity and 1 − specificity across
all observed haemoglobin values; AUC is obtained by trapezoidal integration.
Its confidence interval uses a rescaled bootstrap resampling PSUs with
replacement within strata, 1000 replicates, with the pre-sorted-value
optimisation so that each replicate is a re-weighting rather than a re-sort.

**Decision-curve analysis.** Net benefit of a triage rule at threshold
probability p_t is

  NB = TP_rate − FP_rate · [p_t / (1 − p_t)]

where TP_rate and FP_rate are design-weighted proportions of the whole child
population. The decision here is whether to perform a haemoglobin test, so p_t
is the probability of anaemia at which a worker is indifferent between testing
and not testing; it is the ratio of the cost of one unnecessary test to the
benefit of detecting one anaemic child. Because a point-of-care test is
inexpensive and a missed anaemic child is not, the policy-relevant range is
low. Net benefit is therefore reported across p_t from 0.01 to 0.50 and
compared with the two reference strategies, test-all and test-none.

**Tests-saved / cases-missed frontier.** For each rule, the design-weighted
proportion of children not tested and the proportion of all anaemic children
missed are plotted against each other, with the test-all and test-none
strategies as anchors.

**State-level analysis.** Rule R1 is repeated within each of the 36 states and
union territories, with design-based confidence intervals, and each state is
classified against the Section 8 bar.

**Per-100,000 translation.** For each rule, results are expressed as tests
performed, anaemic children detected and anaemic children missed per 100,000
children aged 6–59 months, to make the trade-off legible to programme managers.

**Missing data.** Complete-case analysis on the inclusion criteria in Section 5.
The number and characteristics of excluded children are reported, and the
weighted anaemia prevalence in the analytic sample is compared with the
published NFHS-5 national figure as an external consistency check.

## 10. Sensitivity analyses

1. Restriction to non-pregnant mothers only, removing any influence of the
   dual maternal threshold.
2. Reference standard varied to haemoglobin < 10.0 g/dL (moderate or severe
   child anaemia), the subgroup with the strongest clinical claim on a scarce
   test.
3. Restriction to children aged 12–59 months, removing the youngest children in
   whom haemoglobin is most physiologically variable.
4. Unadjusted (non-altitude-adjusted) haemoglobin for both mother and child.

## 11. Pre-specified interpretation rules

To prevent result-driven narration, the interpretation is fixed in advance:

- If **any** rule meets the Section 8 bar nationally, the paper reports that
  rule as a candidate for programme use, with its residual miss rate stated
  explicitly, and recommends prospective field validation before adoption.
- If **no** rule meets the bar nationally but some states do, the paper reports
  the state exceptions and examines whether they are explained by lower child
  anaemia prevalence, and states the general principle governing where triage
  can and cannot work.
- If **no** rule meets the bar anywhere, the paper concludes that household
  triage is not a valid substitute for universal child haemoglobin testing, and
  the policy implication is that the Test step must be resourced and counted
  rather than rationed informally.

In all three cases the paper reports the arithmetic reason for the result — the
relationship between target-condition prevalence and attainable negative
predictive value — rather than presenting the finding as a property of Indian
data alone.

## 12. What this study does not claim

1. It makes no causal claim about maternal anaemia and child anaemia.
2. It infers nothing about iron deficiency, or about any other cause of
   anaemia, from haemoglobin.
3. It does not compare NFHS-5 with any other survey round and makes no trend
   claim. The estimand is rule performance within a single survey.
4. It does not treat NFHS-5 capillary point-of-care haemoglobin as equivalent
   to venous laboratory haemoglobin. The consequences of that measurement
   choice for this analysis are addressed in the limitations.
5. It does not evaluate the effectiveness of any treatment.

## 13. Ethics and registration

The analysis uses de-identified, publicly released survey microdata obtained
under an approved DHS Program data request. It involves no contact with human
participants and no identifiable private information. Under the ICMR National
Ethical Guidelines for Biomedical and Health Research Involving Human
Participants (2017), separate institutional ethics review is not required for
secondary analysis of such data. NFHS-5 itself obtained informed consent from
all respondents, including consent for haemoglobin measurement, and was
approved by the ethics review boards of the International Institute for
Population Sciences and ICF.

Prospective registration is not applicable to a secondary analysis of existing
survey data. In its place, this protocol is cryptographically frozen before the
confirmatory analysis is run.

## 14. Reproducibility

| Artefact | Path |
|---|---|
| Protocol freeze record | `governance/PROTOCOL_FREEZE.txt` |
| Survey estimators | `code/svydesign.py` |
| Dataset construction | `code/build_analysis_dataset.py` |
| Confirmatory analysis | `code/run_triage_analysis.py` |
| Figures | `code/make_figures.py` |
| Package verification | `code/verify_package.py` |
| Derived data checksums | `governance/ARTEFACT_SHA256.txt` |

Source microdata are not redistributed; they must be obtained from The DHS
Program. The pipeline runs from the original ZIP archive without manual steps.
