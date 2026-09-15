"""
model_check.py — cross-check of 04_model.R using statsmodels.

The R model uses nlme::lme with a random intercept per cow and an AR(1)
residual correlation structure. statsmodels MixedLM fits the random intercept
but does not offer AR(1) residuals, so the two are not identical models. That
is the point of running it: the fixed-effect estimates should be close, and
where they differ it is the AR(1) term doing something, not a coding error.

What to expect:

  Coefficients   should agree to roughly two decimal places. A large gap in a
                 coefficient means the model formulas have drifted apart.
  Standard errors should be SMALLER here than in R. Ignoring the serial
                 correlation leaves the model believing each weekly sample is
                 more independent than it is, which shrinks the standard
                 errors. In this dataset the R model estimates phi around 0.5,
                 so the difference is not cosmetic.

This second point is the reason the check is worth running at all. It puts a
number on what the correlation structure is buying, instead of taking it on
faith.

Run from the repository root, after the R pipeline:
    python python/model_check.py
"""

import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


def load():
    try:
        milk = pd.read_csv("data/prepared.csv")
        clusters = pd.read_csv("outputs/cow_clusters.csv")[["cow", "cluster"]]
    except FileNotFoundError as e:
        sys.exit(f"Missing {e.filename}. Run the R pipeline first "
                 f"(analysis/01 through 04).")

    df = milk.merge(clusters, on="cow", how="inner")
    df["cluster"] = df["cluster"].astype("category")
    df["diet"] = df["diet"].astype("category")
    df["wc"] = df["week"] - df["week"].median()
    df["wc2"] = df["wc"] ** 2
    return df.sort_values(["cow", "week"]).reset_index(drop=True)


def main():
    df = load()
    print(f"{len(df)} observations, {df['cow'].nunique()} cows, "
          f"{df['cluster'].nunique()} clusters\n")

    model = smf.mixedlm(
        "protein ~ C(cluster) * (wc + wc2) + C(diet)",
        data=df,
        groups=df["cow"],
    )
    fit = model.fit(reml=True)

    print("=== Fixed effects (random intercept only, no AR(1)) ===")
    tbl = pd.DataFrame({
        "coef": fit.params,
        "se": fit.bse,
        "z": fit.tvalues,
        "p": fit.pvalues,
    })
    print(tbl.to_string(float_format=lambda v: f"{v:.4f}"))

    print("\n=== Variance components ===")
    print(f"Cow intercept variance: {float(fit.cov_re.iloc[0, 0]):.5f}")
    print(f"Residual variance:      {fit.scale:.5f}")

    # --- Residual autocorrelation --------------------------------------------
    # The quantity the R model estimates explicitly as phi. Computing it from
    # the residuals here shows whether leaving it unmodelled was defensible.
    df["resid"] = fit.resid
    lag1 = []
    for _, g in df.groupby("cow", observed=True):
        r = g.sort_values("week")["resid"].to_numpy()
        if len(r) > 3:
            r = r - r.mean()
            denom = np.sum(r ** 2)
            if denom > 0:
                lag1.append(np.sum(r[:-1] * r[1:]) / denom)

    lag1 = np.array(lag1)
    print("\n=== Lag-1 autocorrelation of residuals, per cow ===")
    print(f"Mean:   {lag1.mean():.3f}")
    print(f"Median: {np.median(lag1):.3f}")
    print(f"Range:  {lag1.min():.3f} to {lag1.max():.3f}")

    if abs(lag1.mean()) > 0.2:
        print(
            "\nResiduals are clearly autocorrelated. A model without an AR(1)\n"
            "term — this one — reports standard errors that are too small and\n"
            "p-values that are too optimistic. The R model with corAR1 is the\n"
            "one to report; this fit is a cross-check on the coefficients only."
        )
    else:
        print("\nLittle residual autocorrelation; the two models should agree closely.")

    out = tbl.reset_index().rename(columns={"index": "term"})
    out.to_csv("outputs/python_model_coefficients.csv", index=False)
    print("\nWrote outputs/python_model_coefficients.csv")
    print("Compare against outputs/model_summary.txt from the R pipeline.")


if __name__ == "__main__":
    main()
