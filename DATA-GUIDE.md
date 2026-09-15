# Data guide

Everything you need to run this pipeline, understand what each step is
checking, and point it at a different dataset — including a client's.

---

## 1. Install

**R** (4.0 or newer). The pipeline deliberately uses only base R plus packages
that ship with it — `nlme`, `cluster`, `stats` — with `ggplot2` as the single
external dependency. That decision is worth keeping when you adapt this: every
extra package is another thing that fails to install on a slow connection or a
locked-down machine.

```r
install.packages("ggplot2")
```

**Python** (3.9 or newer), only for the cross-check scripts:

```bash
pip install pandas numpy scikit-learn statsmodels
```

---

## 2. Run

From the repository root, in order. Each script writes files the next one
reads, so do not skip.

```bash
Rscript analysis/01_prepare.R
Rscript analysis/02_features.R
Rscript analysis/03_cluster.R
Rscript analysis/04_model.R
Rscript analysis/05_figures.R

python python/cluster_check.py
python python/model_check.py
```

Total runtime is about a minute. The bootstrap in `03` is the slowest part.

---

## 3. What each step does, and what to actually look at

Do not run these and skim past the console output. Each script prints a small
number of things specifically so you can catch a problem before it propagates.

### 01_prepare.R — load and sanity-check

Loads `nlme::Milk` and writes `data/prepared.rds` and `data/prepared.csv`.

**Check:**
- Row count, cow count and week range match what you expect
- Missing values — the column sums should be zero here
- Observations per cow. This dataset is unbalanced, 12 to 19 weeks per cow.
  Unbalanced is fine; knowing it is unbalanced is what matters, because it
  rules out methods that assume a rectangular matrix
- Any implausible values flagged. Nothing is deleted silently

**Setting to know:** `MIN_OBS <- 10` drops cows with too few observations to
support a curve fit. Change it and rerun to see how sensitive your result is.
If the clusters move a lot when you go from 10 to 12, that is information.

### 02_features.R — reduce each series to shape

Writes `data/features.rds` and `outputs/features.csv`.

**Check the correlation matrix it prints.** This is the most useful output in
the script. Any pair above 0.9 means two features are carrying the same
information, which silently doubles the weight of whatever they measure in the
distance metric.

Two real problems surfaced here during development, both caught by that check:

1. Fitting the quadratic on raw week made the linear and quadratic
   coefficients correlate at -0.94. Not because slope and curvature are the
   same thing, but because the intercept sat at week 0, outside the data.
   Centring week at the middle of the observed window fixed it.
2. After that, `net_change` and `slope_linear` correlated at 0.92 — "where it
   ended minus where it started" and "how fast it was changing" are one
   statement twice. `net_change` was dropped from the clustering set and kept
   in the output file for readability.

Neither was visible by eye. Both would have quietly distorted the clusters.

### 03_cluster.R — cluster, then try to break it

Writes `outputs/cow_clusters.csv`, `k_selection.csv`, `cluster_stability.csv`,
`null_benchmark.csv`.

Three separate questions, in order:

**Which k?** Highest average silhouette across k = 2 to 8.

**Is there structure at all?** The null benchmark. k-means optimises
silhouette, so it returns a positive silhouette even on pure noise — the
absolute number is close to meaningless on its own. The script shuffles each
feature independently across cows 200 times, destroying joint structure while
keeping every marginal distribution intact, and reclusters. Compare your
observed value to that distribution. **This is the single most important
output in the repository.**

**Are the individual clusters real?** Bootstrap Jaccard stability, following
Hennig (2007): resample cows with replacement, recluster, assign every
original cow to its nearest bootstrap centroid, measure overlap. Hennig's
reading of the resulting values:

| Mean Jaccard | Meaning |
|---|---|
| below 0.60 | Dissolved. Do not interpret as a cluster |
| 0.60 - 0.75 | A pattern is there, but not a stable one |
| 0.75 - 0.85 | Stable |
| above 0.85 | Highly stable |

### 04_model.R — model, then test properly

Writes `outputs/model_summary.txt` and `model_predictions.csv`.

Two parts, kept deliberately apart, and the separation is the whole point.

**Part A is descriptive.** The mixed-effects model estimates each cluster's
average trajectory with honest uncertainty. The random intercept accounts for
repeated measurements on the same cow; the AR(1) term accounts for adjacent
weeks resembling each other beyond that.

Part A **cannot be used as evidence that the clusters are different**. The
clusters were built from these same trajectories. Asking whether they differ
in trajectory is circular — the algorithm was built to make them differ, and
any p-value from that comparison is meaningless. If you take one thing from
this repository into client work, take this.

**Part B is the actual test.** Diet was deliberately excluded from the feature
set in `02` so it could be used here. If trajectory type is associated with
diet, that is a connection the clustering could not have manufactured.

**Check:** the AR(1) comparison. If the likelihood ratio test says the
correlation structure is not earning its place, drop it rather than keeping it
for appearances.

### 05_figures.R — six figures

Every figure answers a question a reviewer would ask. Figure 1 shows the whole
herd in grey behind each cluster, because a plot of cluster means alone always
looks more convincing than the data deserves.

### The Python cross-checks

Running the same analysis in a second toolchain is the cheapest error check
available. Different defaults, different RNG, different scaling conventions —
if both land on the same partition, it is a property of the data rather than
of a package.

`cluster_check.py` compares partitions using the adjusted Rand index (1.0 is
identical, 0.0 is chance).

`model_check.py` is not fitting the same model: statsmodels MixedLM has no
AR(1) option. That is deliberate. Coefficients should agree; standard errors
should come out **smaller** in Python, because ignoring serial correlation
makes a model believe each weekly sample is more independent than it is. The
script prints the residual autocorrelation so you can see the size of what
you would have been ignoring.

---

## 4. Pointing this at your own data

The pipeline needs long format: one row per subject per time point.

| Column | What it is |
|---|---|
| `cow` | Subject identifier |
| `week` | Time point, numeric |
| `protein` | The measured response |
| `diet` | A subject-level grouping variable held out of the clustering |

Rename your columns to these in `01_prepare.R` and the rest runs unchanged.

**Five things to change deliberately, not by accident:**

1. **`MIN_OBS`** in `01`. Depends on how many time points your curve fit needs.
2. **The curve form** in `02`. A quadratic suits protein percentage over 19
   weeks. For milk yield across a full lactation, a Wilmink or Wood curve is
   the standard choice and will fit the early peak far better.
3. **`CENTRE_WEEK`** in `02`. Should sit inside your observed window.
4. **The feature set** in `02`. The deviation features assume repeated
   measurements dense enough for a wobble to mean something. Monthly test-day
   records will not support `dev_acf1` in the way weekly records do.
5. **The held-out variable** in `04`. You need *something* that was not used
   in clustering, or you have no non-circular test at all. Parity, breed,
   farm, treatment group, season. Decide this before you cluster, not after.

**The one thing not to change:** keep the null benchmark and the stability
test. They are what separate an analysis from a picture.

---

## 5. Common failures

**`contrasts can be applied only to factors with 2 or more levels`**
You built a prediction grid holding a factor at one level. The column still
needs all the levels the model was fitted with:
`factor(levels(dat$diet)[1], levels = levels(dat$diet))`.

**`Error in kmeans: empty cluster`**
k is too large for your sample, or duplicate rows are collapsing. Lower k, or
raise `nstart`.

**Stability comes back near zero for every cluster**
Usually correct rather than broken. It means the variation is continuous and
k-means was partitioning a cloud. Report that. It is a finding.

**Silhouette rises monotonically with k**
Almost always a sign that k-means is splitting noise. Check the null
benchmark before choosing a large k.

**`lme` fails to converge**
Too many fixed effects for the number of subjects, or a random effects
structure the data cannot support. Simplify the random effects before the
fixed ones.

---

## 6. Reading the result this repository produces

The current run gives a result that would be easy to misreport in either
direction, so it is worth stating plainly.

Average silhouette at the chosen k is **0.194**. Against Rousseeuw's rules of
thumb, anything under 0.25 means no substantial structure — so the naive
reading is "nothing here."

But the null benchmark puts shuffled data at a median of **0.132** and a 95th
percentile of **0.145**, giving an empirical p of **0.005**. The observed
separation is well outside what structureless data produces. Bootstrap
stability then runs **0.62 to 0.79**, so the partition reproduces under
resampling.

Both readings are true at once, and the honest summary is: *the trajectory
types are reproducible and not an artefact, but they are regions of a
continuum rather than discrete types.* Cows near a boundary could reasonably
belong to either neighbour.

And the held-out test returns nothing: trajectory type shows no association
with diet (p = 0.79). The clusters describe real, reproducible variation in
trajectory shape, and this analysis provides no evidence that they track diet.

Reporting only the stability numbers would oversell it. Reporting only the
silhouette would bury a real finding. Reporting only Part A would be circular.
Say all three.
