## 01_prepare.R ---------------------------------------------------------------
## Load the public dataset, check its shape, and write a prepared long-format
## table for everything downstream.
##
## Data: nlme::Milk — protein content of milk samples from 79 cows, measured
## weekly from week 1 to week 19 of lactation, across three diet treatments.
## Ships with the nlme package, so this runs with no data request and no
## restricted data in the repository.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(nlme))

dir.create("data", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

data(Milk, package = "nlme")

milk <- data.frame(
  cow     = factor(as.character(Milk$Cow)),   # drop the ordered-factor class
  diet    = factor(as.character(Milk$Diet)),
  week    = as.numeric(Milk$Time),
  protein = as.numeric(Milk$protein)
)

## --- Inspect before trusting -------------------------------------------------
## Everything below is printed rather than assumed. If you swap in your own
## data, this is the block that tells you whether the rest of the pipeline is
## going to behave.

cat("\n=== Structure ===\n")
cat("Rows:               ", nrow(milk), "\n")
cat("Cows:               ", nlevels(milk$cow), "\n")
cat("Week range:         ", paste(range(milk$week), collapse = " to "), "\n")
cat("Diets:              ", paste(levels(milk$diet), collapse = ", "), "\n")

obs_per_cow <- as.numeric(table(milk$cow))
cat("\n=== Balance ===\n")
cat("Observations per cow (min / median / max): ",
    paste(c(min(obs_per_cow), median(obs_per_cow), max(obs_per_cow)),
          collapse = " / "), "\n")
cat("Perfectly balanced:  ", length(unique(obs_per_cow)) == 1, "\n")

cat("\n=== Missing values ===\n")
print(colSums(is.na(milk)))

cat("\n=== Response distribution ===\n")
print(summary(milk$protein))

## --- Quality filters ---------------------------------------------------------
## A cow needs enough points for a trajectory to mean anything. Ten is a
## judgement call, not a rule: with fewer than about ten weekly observations
## the shape features in 02 become noise estimates. Change MIN_OBS and rerun to
## see how sensitive your result is to it.

MIN_OBS <- 10

keep <- names(table(milk$cow))[table(milk$cow) >= MIN_OBS]
dropped <- setdiff(levels(milk$cow), keep)

cat("\n=== Filtering (MIN_OBS =", MIN_OBS, ") ===\n")
cat("Cows kept:    ", length(keep), "\n")
cat("Cows dropped: ", length(dropped),
    if (length(dropped)) paste0(" (", paste(dropped, collapse = ", "), ")") else "", "\n")

milk <- droplevels(milk[milk$cow %in% keep, ])

## Implausible values. For milk protein percentage, anything outside 2-6% is a
## recording error rather than biology. Flag loudly; do not silently delete.
implausible <- milk$protein < 2 | milk$protein > 6
if (any(implausible)) {
  cat("\nWARNING:", sum(implausible), "implausible protein values flagged:\n")
  print(milk[implausible, ])
} else {
  cat("No implausible values.\n")
}

## --- Save --------------------------------------------------------------------
saveRDS(milk, "data/prepared.rds")
write.csv(milk, "data/prepared.csv", row.names = FALSE)

cat("\nWrote data/prepared.rds and data/prepared.csv —",
    nrow(milk), "rows,", nlevels(milk$cow), "cows.\n")
