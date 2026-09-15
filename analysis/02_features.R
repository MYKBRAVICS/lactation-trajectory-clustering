## 02_features.R --------------------------------------------------------------
## Reduce each cow's time series to a small set of interpretable features.
##
## Why features instead of clustering the raw series:
## clustering raw weekly values mostly recovers overall level — high-protein
## cows in one cluster, low-protein cows in another — which tells you something
## you could have got from a single mean. The question is about shape and
## stability, so the features encode shape and stability explicitly.
##
## Two families of feature:
##   1. Curve shape    — where the trajectory starts, how fast it moves, where
##                       it turns, where it ends.
##   2. Deviations     — how far and how erratically the cow departs from its
##                       own smooth curve. Following Poppe et al. (2020), the
##                       variance, autocorrelation and skewness of deviations
##                       from a fitted lactation curve are used as resilience
##                       indicators: a cow that wobbles and takes a long time
##                       to return to its own curve is behaving differently
##                       from one that tracks it closely.
## -----------------------------------------------------------------------------

milk <- readRDS("data/prepared.rds")

## Reference point for the per-cow curve fit. Set to the middle of the
## observed window so that the linear coefficient means "rate of change at
## mid-lactation" rather than "rate of change at week 0", which is outside the
## data and produces collinear coefficients. See the note in cow_features().
CENTRE_WEEK <- median(milk$week)
cat("Centring week at:", CENTRE_WEEK, "\n")

## Sample skewness. Written out rather than pulled from a package so the
## pipeline depends only on base R and packages that ship with it.
skewness <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3) return(NA_real_)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean((x - mean(x))^3) / s^3
}

## Lag-1 autocorrelation of the deviation series.
lag1_acf <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4) return(NA_real_)
  a <- acf(x, lag.max = 1, plot = FALSE, demean = TRUE)$acf
  as.numeric(a[2])
}

cow_features <- function(d) {
  d <- d[order(d$week), ]

  ## Per-cow quadratic in week. With 12-19 points this is stable, and it gives
  ## an analytic turning point. A Wilmink or Wood curve would be the standard
  ## choice for milk yield; for protein percentage over 19 weeks a quadratic
  ## fits the shape without the convergence problems of a nonlinear fit on a
  ## short series.
  ##
  ## Week is centred at CENTRE_WEEK first. On raw week, the linear and
  ## quadratic coefficients come out correlated at about -0.94 — not because
  ## slope and curvature are the same thing, but because the intercept sits at
  ## week 0, outside the data. Centring moves the reference point into the
  ## middle of the observed window: b1 becomes the rate of change at
  ## mid-lactation and is then close to independent of b2. The correlation
  ## check at the bottom of this script is what surfaced this.
  d$wc <- d$week - CENTRE_WEEK

  fit <- lm(protein ~ wc + I(wc^2), data = d)
  cf  <- coef(fit)
  b1  <- unname(cf["wc"])
  b2  <- unname(cf["I(wc^2)"])

  ## Turning point of the fitted parabola, clamped to the observed window.
  ## Unclamped values mean the cow has no turning point inside the lactation
  ## stage observed — worth knowing, so it is recorded separately.
  tp_raw <- if (is.finite(b2) && b2 != 0) CENTRE_WEEK - b1 / (2 * b2) else NA_real_
  tp     <- min(max(tp_raw, min(d$week)), max(d$week))

  dev <- residuals(fit)

  early <- d$protein[d$week <= 4]
  late  <- d$protein[d$week >= max(d$week) - 4]

  data.frame(
    cow            = d$cow[1],
    diet           = d$diet[1],
    n_obs          = nrow(d),
    level_early    = mean(early),                    # starting level
    level_late     = mean(late),                     # settled level
    net_change     = mean(late) - mean(early),       # overall direction
    slope_linear   = b1,                             # initial rate of change
    curvature      = b2,                             # how sharply it turns
    turning_week   = tp,                             # when it turns
    turning_inside = as.numeric(is.finite(tp_raw) &&
                                tp_raw > min(d$week) &&
                                tp_raw < max(d$week)),
    dev_sd         = sd(dev),                        # magnitude of wobble
    dev_acf1       = lag1_acf(dev),                  # persistence of wobble
    dev_skew       = skewness(dev),                  # asymmetry of wobble
    stringsAsFactors = FALSE
  )
}

feat <- do.call(rbind, lapply(split(milk, milk$cow), cow_features))
rownames(feat) <- NULL

## --- Check the features before clustering on them -----------------------------
cat("\n=== Feature summary ===\n")
print(summary(feat[, setdiff(names(feat), c("cow", "diet"))]))

cat("\n=== Missing feature values ===\n")
print(colSums(is.na(feat)))

## Features used for clustering. n_obs and turning_inside are diagnostics, not
## biology, so they stay out. diet stays out deliberately — it is held back as
## an independent variable to test the clusters against in 04.
## net_change is computed and exported, but is NOT clustered on: it correlates
## with slope_linear at 0.92, because "where it ended minus where it started"
## and "how fast it was changing" are the same statement twice. Keeping both
## would silently give overall direction double weight in the distance metric.
## It stays in the output file because it is the most readable summary of the
## three for anyone looking at the cluster table.
CLUSTER_VARS <- c("level_early", "slope_linear", "curvature",
                  "turning_week", "dev_sd", "dev_acf1", "dev_skew")

cat("\n=== Correlation among clustering features ===\n")
cm <- cor(feat[, CLUSTER_VARS], use = "complete.obs")
print(round(cm, 2))

## Near-duplicate features silently double the weight of whatever they measure.
high <- which(abs(cm) > 0.9 & upper.tri(cm), arr.ind = TRUE)
if (nrow(high)) {
  cat("\nWARNING: feature pairs correlated above 0.9 —",
      "they carry the same information twice in the distance metric:\n")
  for (i in seq_len(nrow(high))) {
    cat("  ", CLUSTER_VARS[high[i, 1]], "<->", CLUSTER_VARS[high[i, 2]],
        " r =", round(cm[high[i, 1], high[i, 2]], 2), "\n")
  }
} else {
  cat("\nNo feature pair correlated above 0.9.\n")
}

complete <- complete.cases(feat[, CLUSTER_VARS])
if (any(!complete)) {
  cat("\nDropping", sum(!complete), "cows with incomplete features.\n")
  feat <- feat[complete, ]
}

saveRDS(list(features = feat, cluster_vars = CLUSTER_VARS), "data/features.rds")
write.csv(feat, "outputs/features.csv", row.names = FALSE)

cat("\nWrote data/features.rds and outputs/features.csv —",
    nrow(feat), "cows,", length(CLUSTER_VARS), "clustering features.\n")
