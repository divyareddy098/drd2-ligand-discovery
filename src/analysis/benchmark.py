"""Virtual-screening benchmark metrics.

Calculates ROC-AUC and enrichment factors for ranked molecular
screening results containing known actives and decoys.
"""

from pathlib import Path

import pandas as pd
from sklearn.metrics import roc_auc_score


def enrichment_factor(
    labels: pd.Series,
    scores: pd.Series,
    fraction: float = 0.01,
    lower_is_better: bool = True,
) -> float:
    """Calculate enrichment factor at a specified library fraction."""

    results = pd.DataFrame({
        "active": labels.astype(int),
        "score": scores.astype(float),
    })

    results = results.sort_values(
        "score",
        ascending=lower_is_better,
    )

    n_total = len(results)
    n_actives = int(results["active"].sum())

    if n_total == 0 or n_actives == 0:
        raise ValueError("Dataset must contain compounds and known actives.")

    n_selected = max(1, int(n_total * fraction))

    top = results.head(n_selected)

    actives_selected = int(top["active"].sum())

    active_fraction_library = n_actives / n_total
    active_fraction_top = actives_selected / n_selected

    return active_fraction_top / active_fraction_library


def evaluate_screening(
    input_file: str,
    lower_is_better: bool = True,
):
    """Evaluate a virtual-screening results file."""

    df = pd.read_csv(input_file)

    required = {"compound_id", "active", "docking_score"}

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}"
        )

    df = df.dropna(
        subset=["active", "docking_score"]
    ).copy()

    df["active"] = df["active"].astype(int)
    df["docking_score"] = df["docking_score"].astype(float)

    n_total = len(df)
    n_actives = int(df["active"].sum())
    n_decoys = n_total - n_actives

    # AutoDock/Vina-like docking:
    # more negative score = better prediction.
    #
    # sklearn expects larger prediction values to indicate
    # stronger evidence for the positive class, so invert
    # scores when lower values are better.
    auc_scores = (
        -df["docking_score"]
        if lower_is_better
        else df["docking_score"]
    )

    roc_auc = roc_auc_score(
        df["active"],
        auc_scores,
    )

    ef_1 = enrichment_factor(
        df["active"],
        df["docking_score"],
        fraction=0.01,
        lower_is_better=lower_is_better,
    )

    ef_5 = enrichment_factor(
        df["active"],
        df["docking_score"],
        fraction=0.05,
        lower_is_better=lower_is_better,
    )

    ef_10 = enrichment_factor(
        df["active"],
        df["docking_score"],
        fraction=0.10,
        lower_is_better=lower_is_better,
    )

    print("=" * 55)
    print("VIRTUAL SCREENING BENCHMARK")
    print("=" * 55)

    print(f"Total compounds : {n_total:,}")
    print(f"Known actives   : {n_actives:,}")
    print(f"Decoys          : {n_decoys:,}")

    print()
    print(f"ROC-AUC : {roc_auc:.3f}")
    print(f"EF1%    : {ef_1:.2f}x")
    print(f"EF5%    : {ef_5:.2f}x")
    print(f"EF10%   : {ef_10:.2f}x")

    return {
        "total_compounds": n_total,
        "actives": n_actives,
        "decoys": n_decoys,
        "roc_auc": roc_auc,
        "ef1": ef_1,
        "ef5": ef_5,
        "ef10": ef_10,
    }


if __name__ == "__main__":

    evaluate_screening(
        "results/docking/benchmark_scores.csv"
    )
