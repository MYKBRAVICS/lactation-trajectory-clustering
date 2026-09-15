## 03_cluster.R ---------------------------------------------------------------
## Cluster cows by trajectory shape, choose k, and then test whether the
## clusters are real.
##
## The third step is the one most analyses skip. k-means always returns k
## clusters. Give it uniform noise and it will hand back four tidy groups with
## no complaint. Silhouette width tells you which k separates best, but "best
## available k" is not the same as "there is structure here". Bootstrap
## stability (Hennig, 2007) is what distinguishes them: if a cluster keeps
## reappearing under resampling, it is a feature of the data; if it dissolves,
## it was a partition of continuous variation.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(cluster))

set.seed(42)

fx <- readRDS("data/features.rds")
feat <- fx$features
CLUSTER_VARS <- fx$cluster_vars

## Standardise. Without this, features on larger numeric scales dominate the
## Euclidean distance for no reason other than their units.
X <- scale(as.matrix(feat[, CLUSTER_VARS]))
n <- nrow(X)

K_RANGE  <- 2:8
N_START  <- 100     # random restarts — guards against local minima
N_BOOT   <- 1000    # bootstrap resamples for stability

## --- Choose k ----------------------------------------------------------------

D <- dist(X)

sel <- data.frame(k = K_RANGE, silhouette = NA_real_, wss = NA_real_)
fits <- list()

for (i in seq_along(K_RANGE)) {
  k <- K_RANGE[i]
  km <- kmeans(X, centers = k, nstart = N_START, algorithm = "Hartigan-Wong")
  fits[[as.character(k)]] <- km
  sel$silhouette[i] <- mean(silhouette(km$cluster, D)[, 3])
  sel$wss[i]        <- km$tot.withinss
}

cat("\n=== Choosing k ===\n")
print(sel, row.names = FALSE)

k_best <- sel$k[which.max(sel$silhouette)]
cat("\nBest average silhouette at k =", k_best,
    "(", round(max(sel$silhouette), 3), ")\n")

## Rousseeuw's rough reading of average silhouette width: above 0.7 is strong
## structure, 0.5-0.7 reasonable, 0.25-0.5 weak and possibly artificial, below
## 0.25 means no substantial structure was found.
s <- max(sel$silhouette)
cat("Interpretation: ",
    if (s > 0.7) "strong structure"
    else if (s > 0.5) "reasonable structure"
    else if (s > 0.25) "weak — could be artificial"
    else "no substantial structure", "\n")

km <- fits[[as.character(k_best)]]

## --- Null benchmark -----------------------------------------------------------
## A silhouette of 0.19 sounds low, but low compared to what? k-means run on
## data with no cluster structure at all still produces a positive silhouette,
## because it is optimising exactly that quantity. So run it on structureless
## data of the same size and shape and see what it returns.
##
## The null here preserves the number of cows, the number of features and the
## marginal distribution of each feature, but destroys any joint structure by
## shuffling each feature independently across cows. Anything the real data
## scores above this is structure; anything at or below it is the algorithm
## doing its job on noise.

N_NULL <- 200
null_sil <- numeric(N_NULL)

for (b in seq_len(N_NULL)) {
  Xn <- apply(X, 2, sample)          # shuffle each column independently
  kn <- kmeans(Xn, centers = k_best, nstart = 20, algorithm = "Hartigan-Wong")
  null_sil[b] <- mean(silhouette(kn$cluster, dist(Xn))[, 3])
}

obs_sil <- max(sel$silhouette)
p_emp   <- (sum(null_sil >= obs_sil) + 1) / (N_NULL + 1)

cat("\n=== Null benchmark (", N_NULL, "column-shuffled datasets, k =", k_best, ") ===\n")
cat("Observed silhouette:      ", round(obs_sil, 3), "\n")
cat("Null median:              ", round(median(null_sil), 3), "\n")
cat("Null 95th percentile:     ", round(quantile(null_sil, 0.95), 3), "\n")
cat("Empirical p:              ", round(p_emp, 4), "\n")
cat("Reading: ",
    if (p_emp < 0.05)
      "the separation exceeds what shuffled data produces — the structure is real."
    else
      "the separation is within the range shuffled data produces — treat the clusters as a convenient partition of continuous variation, not as types.",
    "\n")

write.csv(data.frame(observed = obs_sil, null_median = median(null_sil),
                     null_p95 = quantile(null_sil, 0.95), p_empirical = p_emp),
          "outputs/null_benchmark.csv", row.names = FALSE)

## --- Stability ---------------------------------------------------------------
## Hennig's bootstrap scheme: resample cows with replacement, recluster the
## resample, assign every original cow to its nearest bootstrap centroid, and
## measure how much each original cluster overlaps its best bootstrap match
## using the Jaccard coefficient.

nearest_centroid <- function(X, centers) {
  d <- sapply(seq_len(nrow(centers)), function(j)
    colSums((t(X) - centers[j, ])^2))
  max.col(-d, ties.method = "first")
}

jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0) return(0)
  length(intersect(a, b)) / u
}

boot_stability <- function(X, k, original, B) {
  orig_sets <- split(seq_len(nrow(X)), original)
  out <- matrix(NA_real_, nrow = B, ncol = k)

  for (b in seq_len(B)) {
    idx <- sample(seq_len(nrow(X)), nrow(X), replace = TRUE)
    Xb  <- X[idx, , drop = FALSE]

    kb <- tryCatch(
      kmeans(Xb, centers = k, nstart = 10, algorithm = "Hartigan-Wong"),
      error = function(e) NULL)
    if (is.null(kb)) next

    induced  <- nearest_centroid(X, kb$centers)
    boot_sets <- split(seq_len(nrow(X)), induced)

    for (j in seq_len(k)) {
      out[b, j] <- max(vapply(boot_sets,
                              function(s) jaccard(orig_sets[[j]], s),
                              numeric(1)))
    }
  }
  colMeans(out, na.rm = TRUE)
}

cat("\nRunning", N_BOOT, "bootstrap resamples...\n")
stab <- boot_stability(X, k_best, km$cluster, N_BOOT)

verdict <- function(j) {
  if (j >= 0.85) "highly stable"
  else if (j >= 0.75) "stable"
  else if (j >= 0.60) "pattern present, not highly stable"
  else "unstable — do not interpret"
}

stab_tbl <- data.frame(
  cluster = seq_len(k_best),
  n       = as.numeric(table(km$cluster)),
  jaccard = round(stab, 3),
  verdict = vapply(stab, verdict, character(1))
)

cat("\n=== Cluster stability (mean Jaccard over", N_BOOT, "resamples) ===\n")
print(stab_tbl, row.names = FALSE)

retained <- stab_tbl$cluster[stab_tbl$jaccard >= 0.60]
cat("\nClusters passing the 0.60 threshold and carried into the modelling:",
    if (length(retained)) paste(retained, collapse = ", ") else "none", "\n")
if (length(retained) < k_best) {
  cat("Clusters below threshold are reported but not interpreted as types.\n")
}

## --- Describe the clusters ----------------------------------------------------
## Cluster numbers are arbitrary labels. This table is what turns them into
## something a farmer or a reviewer can argue with.

feat$cluster <- factor(km$cluster)

profile <- aggregate(feat[, CLUSTER_VARS], by = list(cluster = feat$cluster), mean)
cat("\n=== Cluster means (original units) ===\n")
print(round(profile[, -1], 3), row.names = FALSE)

cat("\n=== Cluster means (standardised — how far from the herd average) ===\n")
profile_z <- aggregate(as.data.frame(X), by = list(cluster = feat$cluster), mean)
print(round(profile_z[, -1], 2), row.names = FALSE)

saveRDS(list(features = feat, km = km, k = k_best, selection = sel,
             stability = stab_tbl, X = X, retained = retained),
        "data/clusters.rds")

write.csv(feat, "outputs/cow_clusters.csv", row.names = FALSE)
write.csv(sel, "outputs/k_selection.csv", row.names = FALSE)
write.csv(stab_tbl, "outputs/cluster_stability.csv", row.names = FALSE)

cat("\nWrote outputs/cow_clusters.csv, k_selection.csv, cluster_stability.csv\n")
