"""Retrieve and prepare chemical structures for the DRD2 ligand dataset."""

from pathlib import Path
import time

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from chembl_webresource_client.new_client import new_client


INPUT_FILE = Path(
    "data/processed/drd2_compound_affinity.csv"
)

OUTPUT_FILE = Path(
    "data/processed/drd2_compounds_prepared.csv"
)


def retrieve_structures(chembl_ids):
    """Retrieve canonical SMILES from ChEMBL."""

    molecule_client = new_client.molecule

    structures = {}

    total = len(chembl_ids)

    print(
        f"Retrieving structures for {total:,} compounds..."
    )

    batch_size = 100

    for start in range(0, total, batch_size):

        batch = chembl_ids[
            start:start + batch_size
        ]

        try:

            records = molecule_client.filter(
                molecule_chembl_id__in=batch
            ).only(
                [
                    "molecule_chembl_id",
                    "molecule_structures",
                ]
            )

            for record in records:

                chembl_id = record.get(
                    "molecule_chembl_id"
                )

                structures_data = record.get(
                    "molecule_structures"
                )

                if not structures_data:
                    continue

                smiles = structures_data.get(
                    "canonical_smiles"
                )

                if smiles:
                    structures[chembl_id] = smiles

        except Exception as exc:

            print(
                f"Warning: batch starting at "
                f"{start:,} failed: {exc}"
            )

        completed = min(
            start + batch_size,
            total
        )

        print(
            f"\rRetrieved {completed:,}/{total:,}",
            end="",
            flush=True,
        )

        time.sleep(0.05)

    print()

    return structures


def prepare_structure(smiles):
    """Parse and canonicalize a SMILES string using RDKit."""

    if not isinstance(smiles, str):
        return None

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    # Keep largest fragment if multiple disconnected
    # components are present (e.g. salts).
    fragments = Chem.GetMolFrags(
        mol,
        asMols=True,
        sanitizeFrags=True,
    )

    if not fragments:
        return None

    mol = max(
        fragments,
        key=lambda fragment:
            fragment.GetNumHeavyAtoms()
    )

    canonical_smiles = Chem.MolToSmiles(
        mol,
        canonical=True,
    )

    return mol, canonical_smiles


def calculate_descriptors(mol):
    """Calculate basic drug-like molecular descriptors."""

    return {
        "molecular_weight":
            Descriptors.MolWt(mol),

        "logp":
            Descriptors.MolLogP(mol),

        "tpsa":
            Descriptors.TPSA(mol),

        "hbd":
            Descriptors.NumHDonors(mol),

        "hba":
            Descriptors.NumHAcceptors(mol),

        "rotatable_bonds":
            Descriptors.NumRotatableBonds(mol),

        "ring_count":
            Descriptors.RingCount(mol),

        "fraction_csp3":
            Descriptors.FractionCSP3(mol),
    }


def main():

    print(f"Reading affinity dataset: {INPUT_FILE}")

    df = pd.read_csv(INPUT_FILE)

    chembl_ids = (
        df["molecule_chembl_id"]
        .dropna()
        .unique()
        .tolist()
    )

    print(
        f"Unique ChEMBL IDs: {len(chembl_ids):,}"
    )

    structures = retrieve_structures(
        chembl_ids
    )

    print(
        f"Structures retrieved: "
        f"{len(structures):,}"
    )

    rows = []

    invalid_smiles = 0

    for _, row in df.iterrows():

        chembl_id = row[
            "molecule_chembl_id"
        ]

        smiles = structures.get(
            chembl_id
        )

        if smiles is None:
            continue

        prepared = prepare_structure(
            smiles
        )

        if prepared is None:
            invalid_smiles += 1
            continue

        mol, canonical_smiles = prepared

        descriptors = calculate_descriptors(
            mol
        )

        output = row.to_dict()

        output[
            "original_smiles"
        ] = smiles

        output[
            "canonical_smiles"
        ] = canonical_smiles

        output.update(descriptors)

        rows.append(output)

    prepared_df = pd.DataFrame(rows)

    before_duplicates = len(
        prepared_df
    )

    prepared_df = (
        prepared_df
        .sort_values(
            "median_ki_nm"
        )
        .drop_duplicates(
            subset="canonical_smiles",
            keep="first",
        )
        .reset_index(drop=True)
    )

    duplicate_structures = (
        before_duplicates
        - len(prepared_df)
    )

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    prepared_df.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    print("\nCHEMICAL STRUCTURE QC")
    print("=" * 50)

    print(
        f"Starting compounds: "
        f"{len(df):,}"
    )

    print(
        f"Structures retrieved: "
        f"{len(structures):,}"
    )

    print(
        f"Invalid SMILES: "
        f"{invalid_smiles:,}"
    )

    print(
        f"Duplicate standardized structures: "
        f"{duplicate_structures:,}"
    )

    print(
        f"Final unique structures: "
        f"{len(prepared_df):,}"
    )

    print("\nFinal potency classes")
    print("-" * 50)

    print(
        "Median Ki <= 10 nM:",
        f"{(prepared_df['median_ki_nm'] <= 10).sum():,}"
    )

    print(
        "Median Ki <= 100 nM:",
        f"{(prepared_df['median_ki_nm'] <= 100).sum():,}"
    )

    print(
        "Median Ki <= 1,000 nM:",
        f"{(prepared_df['median_ki_nm'] <= 1000).sum():,}"
    )

    print(
        "Median Ki > 1,000 nM:",
        f"{(prepared_df['median_ki_nm'] > 1000).sum():,}"
    )

    print(
        f"\nSaved prepared dataset to: "
        f"{OUTPUT_FILE}"
    )


if __name__ == "__main__":
    main()
