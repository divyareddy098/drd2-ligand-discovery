"""Build a property-matched DRD2 weak-binder benchmark set."""

from pathlib import Path

import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


COMPOUND_FILE = Path(
    "data/processed/drd2_compounds_prepared.csv"
)

ACTIVE_FILE = Path(
    "data/processed/drd2_benchmark_actives.csv"
)

OUTPUT_FILE = Path(
    "data/processed/drd2_benchmark_decoys.csv"
)

N_DECOYS = 5000

WEAK_KI_THRESHOLD = 1000.0

MAX_ACTIVE_SIMILARITY = 0.50

RADIUS = 2
FP_SIZE = 2048

RANDOM_SEED = 42


def fingerprint_generator():

    return rdFingerprintGenerator.GetMorganGenerator(
        radius=RADIUS,
        fpSize=FP_SIZE,
    )


def make_fingerprints(smiles, generator):

    molecules = [
        Chem.MolFromSmiles(s)
        for s in smiles
    ]

    return [
        generator.GetFingerprint(m)
        for m in molecules
    ]


def property_filter(df, actives):

    """Restrict candidates to the broad property space of actives."""

    properties = [
        "molecular_weight",
        "logp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
    ]

    filtered = df.copy()

    print("\nActive property ranges")
    print("-" * 50)

    for prop in properties:

        low = actives[prop].quantile(0.01)
        high = actives[prop].quantile(0.99)

        print(
            f"{prop:18s}: "
            f"{low:.2f} - {high:.2f}"
        )

        filtered = filtered[
            (filtered[prop] >= low)
            & (filtered[prop] <= high)
        ]

    return filtered


def remove_active_like_candidates(
    candidates,
    active_fps,
    generator,
):

    """Remove weak binders highly similar to benchmark actives."""

    kept_rows = []
    maximum_similarities = []

    total = len(candidates)

    print(
        "\nRemoving candidates structurally similar "
        "to actives..."
    )

    for count, (_, row) in enumerate(
        candidates.iterrows(),
        start=1,
    ):

        mol = Chem.MolFromSmiles(
            row["canonical_smiles"]
        )

        if mol is None:
            continue

        fp = generator.GetFingerprint(mol)

        similarities = (
            DataStructs.BulkTanimotoSimilarity(
                fp,
                active_fps,
            )
        )

        max_similarity = max(
            similarities
        )

        if max_similarity < MAX_ACTIVE_SIMILARITY:

            kept_rows.append(
                row
            )

            maximum_similarities.append(
                max_similarity
            )

        if count % 250 == 0:

            print(
                f"\rProcessed "
                f"{count:,}/{total:,}",
                end="",
                flush=True,
            )

    print()

    result = pd.DataFrame(
        kept_rows
    ).reset_index(drop=True)

    result[
        "max_active_tanimoto"
    ] = maximum_similarities

    return result


def main():

    compounds = pd.read_csv(
        COMPOUND_FILE
    )

    actives = pd.read_csv(
        ACTIVE_FILE
    )

    print("=" * 60)
    print("DRD2 BENCHMARK NEGATIVE/DECOY CONSTRUCTION")
    print("=" * 60)

    print(
        f"Prepared structures: "
        f"{len(compounds):,}"
    )

    print(
        f"Benchmark actives: "
        f"{len(actives):,}"
    )

    # Experimentally weak DRD2 binders
    candidates = compounds[
        compounds["median_ki_nm"]
        > WEAK_KI_THRESHOLD
    ].copy()

    print(
        f"Weak-binder candidates "
        f"(Ki > {WEAK_KI_THRESHOLD:.0f} nM): "
        f"{len(candidates):,}"
    )

    # Remove any active benchmark structures
    active_smiles = set(
        actives["canonical_smiles"]
    )

    candidates = candidates[
        ~candidates[
            "canonical_smiles"
        ].isin(active_smiles)
    ].copy()

    candidates = property_filter(
        candidates,
        actives,
    )

    print(
        f"\nAfter property matching: "
        f"{len(candidates):,}"
    )

    generator = fingerprint_generator()

    print(
        "\nGenerating active fingerprints..."
    )

    active_fps = make_fingerprints(
        actives["canonical_smiles"].tolist(),
        generator,
    )

    candidates = remove_active_like_candidates(
        candidates,
        active_fps,
        generator,
    )

    print(
        f"After similarity exclusion "
        f"(Tanimoto < {MAX_ACTIVE_SIMILARITY}): "
        f"{len(candidates):,}"
    )

    # Randomize after all scientific filters.
    candidates = candidates.sample(
        frac=1,
        random_state=RANDOM_SEED,
    ).reset_index(drop=True)

    if len(candidates) > N_DECOYS:

        decoys = candidates.head(
            N_DECOYS
        ).copy()

    else:

        decoys = candidates.copy()

    decoys[
        "benchmark_label"
    ] = 0

    decoys[
        "benchmark_class"
    ] = "weak_binder"

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    decoys.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    print("\nFINAL BENCHMARK NEGATIVE SET")
    print("=" * 60)

    print(
        f"Selected weak-binder controls: "
        f"{len(decoys):,}"
    )

    print(
        f"Median Ki: "
        f"{decoys['median_ki_nm'].median():,.2f} nM"
    )

    print(
        f"Median MW: "
        f"{decoys['molecular_weight'].median():.1f}"
    )

    print(
        f"Median LogP: "
        f"{decoys['logp'].median():.2f}"
    )

    print(
        f"Median max-active Tanimoto: "
        f"{decoys['max_active_tanimoto'].median():.3f}"
    )

    print(
        f"Maximum active similarity: "
        f"{decoys['max_active_tanimoto'].max():.3f}"
    )

    print(
        f"\nSaved: {OUTPUT_FILE}"
    )


if __name__ == "__main__":
    main()
