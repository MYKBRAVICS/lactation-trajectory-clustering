## 04_model.R -----------------------------------------------------------------
## Two separate jobs, deliberately kept apart.
##
## PART A (descriptive). Estimate the average trajectory of each cluster with
## honest uncertainty, using a mixed-effects model that respects the fact that
## the 1337 observations are 79 cows measured repeatedly, not 1337 independent
## measurements.
##
## PART B (inferential). Test the clusters against something that was never
## used to build them.
##
## The distinction matters and is the most common way this kind of analysis
## goes wrong. Clusters here were derived FROM the protein trajectories. Asking
## "do the clusters differ in protein trajectory?" therefore cannot fail — the
## algorithm was built to make them differ, and a p-value from that comparison
## is meaningless. It is circular.
##
## Diet, by contrast, was held out of the clustering entirely. If trajectory
## types distribute unevenly across diets, that is an association the
## clustering could not have manufactured. That is the real test.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(nlme))

set.seed(42)

cl   <- readRDS("data/clusters.rds")
milk <- readRDS("data/prepared.rds")
feat <- cl$features

dat <- merge(milk, feat[, c("cow", "cluster")], by = "cow")
dat$cluster <- factor(dat$cluster)
dat$wc      <- dat$week - median(dat$week)   # same centring as the features
dat         <- dat[order(dat$cow, dat$week), ]

cat("Modelling", nrow(dat), "observations from", nlevels(dat$cow), "cows in",
    nlevels(dat$cluster), "clusters.\n")

## =============================================================================
## PART A — descriptive trajectory model
## =============================================================================

## Random intercept for cow, because repeated measurements on one cow are
## correlated. Ignoring this treats each weekly sample as independent evidence
## and shrinks the standard errors to roughly a quarter of their honest size.
##
## corAR1 on top of that, because adjacent weeks are more alike than distant
## weeks even after accounting for the cow's own level. form = ~ week | cow
## tells nlme the position of each observation in time within each cow, which
## matters here because the series are unbalanced (12-19 weeks per cow).

m_full <- lme(protein ~ cluster * (wc + I(wc^2)) + diet,
              random      = ~ 1 | cow,
              correlation = corAR1(form = ~ week | cow),
              data        = dat,
              method      = "REML")

cat("\n=== Fixed effects ===\n")
print(round(summary(m_full)$tTable, 4))

cat("\n=== Variance components ===\n")
vc <- VarCorr(m_full)
print(vc)
cat("Residual autocorrelation (phi): ",
    round(coef(m_full$modelStruct$corStruct, unconstrained = FALSE), 3), "\n")

## Is the AR1 structure earning its place? Compare against the same model
## without it. If it is not, say so and drop it.
m_noar <- lme(protein ~ cluster * (wc + I(wc^2)) + diet,
              random = ~ 1 | cow, data = dat, method = "REML")

cat("\n=== Does the AR1 correlation structure improve the fit? ===\n")
print(anova(m_noar, m_full))

## --- Model-implied cluster trajectories ---------------------------------------
grid <- expand.grid(
  wc      = seq(min(dat$wc), max(dat$wc), length.out = 60),
  cluster = levels(dat$cluster),
  stringsAsFactors = FALSE
)
grid$cluster <- factor(grid$cluster, levels = levels(dat$cluster))
## Diet is held at its reference level, but the column must carry ALL the
## levels the model was fitted with — otherwise model.matrix cannot build the
## contrasts and fails with "contrasts can be applied only to factors with 2
## or more levels".
grid$diet <- factor(levels(dat$diet)[1], levels = levels(dat$diet))
grid$week <- grid$wc + median(dat$week)

X  <- model.matrix(~ cluster * (wc + I(wc^2)) + diet, data = grid)
b  <- fixef(m_full)
V  <- vcov(m_full)

grid$fit <- as.numeric(X %*% b)
grid$se  <- sqrt(rowSums((X %*% V) * X))
grid$lo  <- grid$fit - 1.96 * grid$se
grid$hi  <- grid$fit + 1.96 * grid$se

write.csv(grid, "outputs/model_predictions.csv", row.names = FALSE)

cat("\n=== Model-implied protein at weeks 2, 10 and 18, by cluster ===\n")
show <- do.call(rbind, lapply(c(2, 10, 18), function(w) {
  g <- grid[abs(grid$week - w) == min(abs(grid$week - w)), ]
  g <- g[!duplicated(g$cluster), ]
  data.frame(week = w, cluster = g$cluster,
             fit = round(g$fit, 3),
             ci = paste0("[", round(g$lo, 3), ", ", round(g$hi, 3), "]"))
}))
print(show, row.names = FALSE)

cat("\nNOTE: the differences above are DESCRIPTIVE. The clusters were built\n",
    "from these trajectories, so separation here is guaranteed by construction\n",
    "and carries no inferential weight. Part B is the test.\n")

## =============================================================================
## PART B — the non-circular test
## =============================================================================
## Diet was excluded from the feature set in 02 precisely so it could be used
## here. If trajectory type is associated with diet, the clustering has found
## something that connects to the world outside its own input.

tab <- table(cluster = feat$cluster, diet = feat$diet)

cat("\n=== Cluster by diet ===\n")
print(tab)
cat("\nRow proportions:\n")
print(round(prop.table(tab, 1), 3))

## Cell counts are small, so use a Monte Carlo Fisher test rather than relying
## on the chi-square approximation.
ft <- fisher.test(tab, simulate.p.value = TRUE, B = 20000)
cat("\nFisher exact test (Monte Carlo,", 20000, "replicates)\n")
cat("p =", round(ft$p.value, 4), "\n")

cs <- suppressWarnings(chisq.test(tab))
cat("Chi-square p =", round(cs$p.value, 4),
    "(approximation — check expected counts below)\n")
cat("Minimum expected count:", round(min(cs$expected), 2),
    if (min(cs$expected) < 5) " — below 5, trust the Fisher result\n" else "\n")

cat("\nVERDICT: ",
    if (ft$p.value < 0.05)
      "trajectory type is associated with diet. The clusters connect to a variable they were not built from."
    else
      "no detectable association between trajectory type and diet. The clusters describe variation in trajectory shape, but nothing here shows they track an external factor.",
    "\n")

## --- Save ---------------------------------------------------------------------
sink("outputs/model_summary.txt")
cat("PART A — trajectory model (descriptive)\n\n")
print(summary(m_full))
cat("\n\nPART B — cluster by diet (non-circular test)\n\n")
print(tab)
cat("\nFisher exact (Monte Carlo) p =", ft$p.value, "\n")
sink()

saveRDS(list(model = m_full, grid = grid, diet_table = tab, fisher = ft),
        "data/model.rds")

cat("\nWrote outputs/model_summary.txt and outputs/model_predictions.csv\n")
