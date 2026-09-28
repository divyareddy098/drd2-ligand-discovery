"""Curate DRD2 binding-affinity data from ChEMBL."""

from pathlib import Path

import pandas as pd
from chembl_webresource_client.new_client import new_client


TARGET_ID = "CHEMBL217"

OUTPUT_DIR = Path("data/processed")
OUTPUT_FILE = OUTPUT_DIR / "drd2_ki_measurements.csv"


def retrieve_drd2_ki():
    """Retrieve human DRD2 Ki binding measurements from ChEMBL."""

    print(f"Retrieving DRD2 data from ChEMBL ({TARGET_ID})...")

    activity = new_client.activity

    records = list(
        activity.filter(
            target_chembl_id=TARGET_ID,
            assay_type="B",
            standard_type="Ki",
            standard_units="nM",
        )
    )

    print(f"Retrieved records: {len(records):,}")

    return records


def curate_records(records):
    """Remove records with missing, invalid, or non-positive Ki values."""

    rows = []

    for record in records:

        molecule_id = record.get("molecule_chembl_id")
        value = record.get("standard_value")

        if molecule_id is None or value is None:
            continue

        try:
            value = float(value)
        except (TypeError, ValueError):
            continue

        if value <= 0:
            continue

        rows.append(
            {
                "molecule_chembl_id": molecule_id,
                "ki_nm": value,
                "pchembl_value": record.get("pchembl_value"),
                "assay_chembl_id": record.get("assay_chembl_id"),
                "document_chembl_id": record.get("document_chembl_id"),
            }
        )

    return pd.DataFrame(rows)


def summarize_dataset(df):
    """Print summary statistics for the curated dataset."""

    print("\nDRD2 Ki DATASET")
    print("=" * 45)

    print(f"Valid Ki measurements: {len(df):,}")

    print(
        f"Unique compounds: "
        f"{df['molecule_chembl_id'].nunique():,}"
    )

    print("\nPotency distribution")
    print("-" * 45)

    print(
        f"Ki <= 10 nM: "
        f"{(df['ki_nm'] <= 10).sum():,}"
    )

    print(
        f"Ki <= 100 nM: "
        f"{(df['ki_nm'] <= 100).sum():,}"
    )

    print(
        f"Ki <= 1,000 nM: "
        f"{(df['ki_nm'] <= 1000).sum():,}"
    )

    print(
        f"Ki > 1,000 nM: "
        f"{(df['ki_nm'] > 1000).sum():,}"
    )

    print("\nKi summary (nM)")
    print("-" * 45)

    print(
        df["ki_nm"]
        .describe(
            percentiles=[0.25, 0.50, 0.75, 0.90]
        )
        .round(2)
    )


def main():
    """Run the complete DRD2 Ki data-curation workflow."""

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    records = retrieve_drd2_ki()

    df = curate_records(records)

    df.to_csv(
        OUTPUT_FILE,
        index=False
    )

    summarize_dataset(df)

    print(
        f"\nSaved curated dataset to: {OUTPUT_FILE}"
    )


if __name__ == "__main__":
    main()
