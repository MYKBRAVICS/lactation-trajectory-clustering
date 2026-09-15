"""
cluster_check.py — independent reimplementation of 03_cluster.R in Python.

Why bother reimplementing something that already works in R:

Running the same analysis twice in two toolchains is the cheapest error check
available. R and Python use different k-means defaults, different random
number generators and different scaling conventions. If both land on the same
partition, the result is a property of the data. If they disagree, one of them
encodes an assumption nobody wrote down — and it is better to find that out
here than in review.

Agreement is measured with the adjusted Rand index, which compares two
partitions while correcting for the agreement expected by chance. 1.0 is
identical; 0.0 is what two unrelated partitions would score.

Run from the repository root, after the R pipeline:
    python python/cluster_check.py
"""

import sys

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.preprocessing import StandardScaler

CLUSTER_VARS = [
    "level_early", "slope_linear", "curvature",
    "turning_week", "dev_sd", "dev_acf1", "dev_skew",
]
K_RANGE = range(2, 9)
SEED = 42


def load():
    try:
        feat = pd.read_csv("outputs/features.csv")
        r_assign = pd.read_csv("outputs/cow_clusters.csv")[["cow", "cluster"]]
    except FileNotFoundError as e:
        sys.exit(f"Missing {e.filename}. Run the R pipeline first "
                 f"(analysis/01 through 03).")
    return feat.merge(r_assign, on="cow", how="left")


def main():
    df = load()
    X = StandardScaler().fit_transform(df[CLUSTER_VARS].to_numpy())

    print(f"{len(df)} cows, {len(CLUSTER_VARS)} features\n")
    print("Note: inertia here will not match the R script's WSS exactly. R's\n"
          "scale() divides by the sample SD (n-1 denominator); sklearn's\n"
          "StandardScaler divides by the population SD (n). With n=79 that is\n"
          "about a 1.3% difference in scale, so the sums of squares differ\n"
          "slightly while the partition itself is unaffected. Compare the\n"
          "silhouette values and the Rand index, not the inertia.\n")

    # --- Choose k, same criterion as the R script ----------------------------
    rows = []
    fits = {}
    for k in K_RANGE:
        km = KMeans(n_clusters=k, n_init=100, random_state=SEED).fit(X)
        fits[k] = km
        rows.append({
            "k": k,
            "silhouette": silhouette_score(X, km.labels_),
            "inertia": km.inertia_,
        })
    sel = pd.DataFrame(rows)
    print("=== Choosing k ===")
    print(sel.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    k_best = int(sel.loc[sel["silhouette"].idxmax(), "k"])
    print(f"\nPython best k: {k_best} "
          f"(silhouette {sel['silhouette'].max():.3f})")

    py_labels = fits[k_best].labels_
    r_labels = df["cluster"].to_numpy()

    # --- Compare to R ---------------------------------------------------------
    print("\n=== Agreement with the R pipeline ===")
    r_k = int(pd.Series(r_labels).nunique())
    print(f"R chose k = {r_k}, Python chose k = {k_best}")

    if r_k != k_best:
        print("MISMATCH on k. The partitions below are still compared, but the "
              "disagreement itself is the finding — check nstart/n_init and "
              "the scaling convention in both scripts before trusting either.")

    ari = adjusted_rand_score(r_labels, py_labels)
    print(f"Adjusted Rand index: {ari:.3f}")

    if ari > 0.95:
        verdict = "identical partition (label numbering aside)"
    elif ari > 0.80:
        verdict = "substantially the same partition; a handful of boundary cows differ"
    elif ari > 0.50:
        verdict = "broadly similar but materially different — investigate before reporting"
    else:
        verdict = "the two toolchains disagree. Do not report either until you know why"
    print(f"Reading: {verdict}")

    # Cross-tabulation makes the disagreement concrete rather than a single number.
    print("\nR clusters (rows) against Python clusters (columns):")
    print(pd.crosstab(pd.Series(r_labels, name="R"),
                      pd.Series(py_labels + 1, name="Python")))

    # --- Null benchmark, same logic as the R script ---------------------------
    rng = np.random.default_rng(SEED)
    null = []
    for _ in range(200):
        Xn = np.column_stack([rng.permutation(col) for col in X.T])
        kn = KMeans(n_clusters=k_best, n_init=20, random_state=SEED).fit(Xn)
        null.append(silhouette_score(Xn, kn.labels_))
    null = np.array(null)

    obs = sel["silhouette"].max()
    p_emp = (np.sum(null >= obs) + 1) / (len(null) + 1)
    print("\n=== Null benchmark (200 column-shuffled datasets) ===")
    print(f"Observed silhouette:  {obs:.3f}")
    print(f"Null median:          {np.median(null):.3f}")
    print(f"Null 95th percentile: {np.quantile(null, 0.95):.3f}")
    print(f"Empirical p:          {p_emp:.4f}")

    out = df[["cow", "cluster"]].copy()
    out["python_cluster"] = py_labels + 1
    out.to_csv("outputs/cluster_agreement.csv", index=False)
    print("\nWrote outputs/cluster_agreement.csv")


if __name__ == "__main__":
    main()
