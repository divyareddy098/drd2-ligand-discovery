"""Aggregate repeated DRD2 Ki measurements to compound-level values."""

from pathlib import Path

import pandas as pd


INPUT_FILE = Path(
    "data/processed/drd2_ki_measurements.csv"
)

OUTPUT_FILE = Path(
    "data/processed/drd2_compound_affinity.csv"
)


def aggregate_compounds(df):
    """Calculate compound-level DRD2 affinity statistics."""

    compound_df = (
        df.groupby("molecule_chembl_id")
        .agg(
            median_ki_nm=("ki_nm", "median"),
            mean_ki_nm=("ki_nm", "mean"),
            min_ki_nm=("ki_nm", "min"),
            max_ki_nm=("ki_nm", "max"),
            n_measurements=("ki_nm", "count"),
        )
        .reset_index()
    )

    return compound_df


def summarize(df):

    print("\nCOMPOUND-LEVEL DRD2 DATASET")
    print("=" * 50)

    print(f"Unique compounds: {len(df):,}")

    print("\nMedian Ki potency classes")
    print("-" * 50)

    very_high = (
        df["median_ki_nm"] <= 10
    ).sum()

    high = (
        df["median_ki_nm"] <= 100
    ).sum()

    moderate = (
        df["median_ki_nm"] <= 1000
    ).sum()

    weak = (
        df["median_ki_nm"] > 1000
    ).sum()

    print(f"Median Ki <= 10 nM:    {very_high:,}")
    print(f"Median Ki <= 100 nM:   {high:,}")
    print(f"Median Ki <= 1,000 nM: {moderate:,}")
    print(f"Median Ki > 1,000 nM:  {weak:,}")

    print("\nMeasurement replication")
    print("-" * 50)

    print(
        "Compounds with >=2 measurements:",
        f"{(df['n_measurements'] >= 2).sum():,}"
    )

    print(
        "Compounds with >=3 measurements:",
        f"{(df['n_measurements'] >= 3).sum():,}"
    )

    print(
        "Maximum measurements for one compound:",
        df["n_measurements"].max()
    )

    print("\nCompound-level median Ki summary")
    print("-" * 50)

    print(
        df["median_ki_nm"]
        .describe(
            percentiles=[0.25, 0.50, 0.75, 0.90]
        )
        .round(2)
    )


def main():

    print(f"Reading: {INPUT_FILE}")

    df = pd.read_csv(INPUT_FILE)

    # Explicit QC
    df = df[
        df["ki_nm"].notna()
        & (df["ki_nm"] > 0)
    ].copy()

    print(
        f"Valid positive measurements: {len(df):,}"
    )

    compounds = aggregate_compounds(df)

    compounds.to_csv(
        OUTPUT_FILE,
        index=False
    )

    summarize(compounds)

    print(
        f"\nSaved compound dataset to: {OUTPUT_FILE}"
    )


if __name__ == "__main__":
    main()
