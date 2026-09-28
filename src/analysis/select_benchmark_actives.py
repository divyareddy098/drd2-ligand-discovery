"""Select a chemically diverse DRD2 active set for docking benchmark."""

from pathlib import Path

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs


INPUT_FILE = Path(
    "data/processed/drd2_compounds_prepared.csv"
)

OUTPUT_FILE = Path(
    "data/processed/drd2_benchmark_actives.csv"
)

KI_THRESHOLD = 100.0
N_ACTIVES = 500

# Morgan fingerprint parameters
RADIUS = 2
FP_SIZE = 2048


def prepare_active_pool(df):
    """Select potent, drug-like, RDKit-valid DRD2 compounds."""

    print(
        f"Starting compounds: {len(df):,}"
    )

    df = df[
        df["median_ki_nm"] <= KI_THRESHOLD
    ].copy()

    print(
        f"Ki <= {KI_THRESHOLD:.0f} nM: "
        f"{len(df):,}"
    )

    df = df[
        (df["molecular_weight"] >= 150)
        & (df["molecular_weight"] <= 600)
        & (df["logp"] >= -2)
        & (df["logp"] <= 6)
        & (df["tpsa"] <= 150)
    ].copy()

    print(
        f"After property QC: {len(df):,}"
    )

    molecules = []
    valid_rows = []

    for index, row in df.iterrows():

        mol = Chem.MolFromSmiles(
            row["canonical_smiles"]
        )

        if mol is not None:

            molecules.append(mol)
            valid_rows.append(index)

    df = df.loc[
        valid_rows
    ].reset_index(drop=True)

    print(
        f"RDKit-valid structures: {len(df):,}"
    )

    return df, molecules


def calculate_fingerprints(molecules):
    """Generate Morgan fingerprints."""

    generator = (
        rdFingerprintGenerator.GetMorganGenerator(
            radius=RADIUS,
            fpSize=FP_SIZE,
        )
    )

    fingerprints = [
        generator.GetFingerprint(mol)
        for mol in molecules
    ]

    return fingerprints


def diverse_selection(
    df,
    fingerprints,
    n_select,
):
    """Select compounds using greedy max-min chemical diversity."""

    if n_select > len(df):
        n_select = len(df)

    # Start with the most potent compound.
    first_index = (
        df["median_ki_nm"]
        .astype(float)
        .idxmin()
    )

    selected = [first_index]

    remaining = set(
        range(len(df))
    )

    remaining.remove(
        first_index
    )

    # For every candidate, keep track of its
    # maximum similarity to anything already selected.
    max_similarity = np.zeros(
        len(df),
        dtype=float,
    )

    print(
        f"\nSelecting {n_select:,} "
        f"chemically diverse actives..."
    )

    while len(selected) < n_select:

        newest = selected[-1]

        candidates = list(
            remaining
        )

        similarities = (
            DataStructs.BulkTanimotoSimilarity(
                fingerprints[newest],
                [
                    fingerprints[i]
                    for i in candidates
                ],
            )
        )

        for candidate, similarity in zip(
            candidates,
            similarities,
        ):
            if similarity > max_similarity[candidate]:
                max_similarity[candidate] = similarity

        # Choose the molecule least similar to
        # anything already selected.
        next_index = min(
            remaining,
            key=lambda i: max_similarity[i],
        )

        selected.append(
            next_index
        )

        remaining.remove(
            next_index
        )

        if (
            len(selected) % 50 == 0
            or len(selected) == n_select
        ):

            print(
                f"Selected "
                f"{len(selected):,}/{n_select:,}"
            )

    selected_df = (
        df.iloc[selected]
        .copy()
        .reset_index(drop=True)
    )

    return selected_df


def similarity_summary(
    selected_df,
    fingerprints,
    original_df,
):
    """Calculate pairwise similarity statistics."""

    selected_smiles = set(
        selected_df["canonical_smiles"]
    )

    selected_indices = [
        i
        for i, smiles in enumerate(
            original_df["canonical_smiles"]
        )
        if smiles in selected_smiles
    ]

    pairwise = []

    for position, i in enumerate(
        selected_indices
    ):

        if position + 1 >= len(
            selected_indices
        ):
            break

        others = selected_indices[
            position + 1:
        ]

        similarities = (
            DataStructs.BulkTanimotoSimilarity(
                fingerprints[i],
                [
                    fingerprints[j]
                    for j in others
                ],
            )
        )

        pairwise.extend(
            similarities
        )

    pairwise = np.asarray(
        pairwise,
        dtype=float,
    )

    return {
        "median": np.median(pairwise),
        "mean": np.mean(pairwise),
        "p90": np.percentile(pairwise, 90),
        "maximum": np.max(pairwise),
    }


def main():

    df = pd.read_csv(
        INPUT_FILE
    )

    active_pool, molecules = (
        prepare_active_pool(df)
    )

    print(
        "\nGenerating Morgan fingerprints..."
    )

    fingerprints = (
        calculate_fingerprints(
            molecules
        )
    )

    print(
        f"Fingerprints generated: "
        f"{len(fingerprints):,}"
    )

    benchmark = diverse_selection(
        active_pool,
        fingerprints,
        N_ACTIVES,
    )

    benchmark[
        "benchmark_label"
    ] = 1

    benchmark.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    stats = similarity_summary(
        benchmark,
        fingerprints,
        active_pool,
    )

    print("\nFINAL BENCHMARK ACTIVE SET")
    print("=" * 50)

    print(
        f"Selected compounds: "
        f"{len(benchmark):,}"
    )

    print(
        f"Median Ki: "
        f"{benchmark['median_ki_nm'].median():.2f} nM"
    )

    print(
        f"Ki range: "
        f"{benchmark['median_ki_nm'].min():.4g}"
        f" - "
        f"{benchmark['median_ki_nm'].max():.2f} nM"
    )

    print("\nPairwise Tanimoto similarity")
    print("-" * 50)

    print(
        f"Mean:    "
        f"{stats['mean']:.3f}"
    )

    print(
        f"Median:  "
        f"{stats['median']:.3f}"
    )

    print(
        f"90th %:  "
        f"{stats['p90']:.3f}"
    )

    print(
        f"Maximum: "
        f"{stats['maximum']:.3f}"
    )

    print(
        f"\nSaved: {OUTPUT_FILE}"
    )


if __name__ == "__main__":
    main()
