"""Repair missing heavy atoms in the chain-break-aware 6CM4 DRD2 receptor."""

from pathlib import Path

from pdbfixer import PDBFixer
from openmm.app import PDBFile


INPUT = Path(
    "structures/processed/6CM4_DRD2_split.pdb"
)

OUTPUT = Path(
    "structures/processed/6CM4_DRD2_split_repaired.pdb"
)


def main():

    print("=" * 60)
    print("6CM4 SPLIT DRD2 RECEPTOR REPAIR")
    print("=" * 60)

    print(f"Input: {INPUT}")

    fixer = PDBFixer(
        filename=str(INPUT)
    )

    fixer.findNonstandardResidues()

    print(
        f"Nonstandard residues: "
        f"{len(fixer.nonstandardResidues)}"
    )

    # Initialize missing-residue information.
    fixer.findMissingResidues()

    detected_missing = dict(
        fixer.missingResidues
    )

    print(
        f"Missing residue segments detected: "
        f"{len(detected_missing)}"
    )

    # Do not reconstruct completely unresolved residues.
    fixer.missingResidues = {}

    # Find missing atoms in residues that are actually present.
    fixer.findMissingAtoms()

    n_missing_atoms = sum(
        len(atoms)
        for atoms in fixer.missingAtoms.values()
    )

    n_terminal_atoms = sum(
        len(atoms)
        for atoms in fixer.missingTerminals.values()
    )

    print()
    print("Atom repair")
    print(
        f"Missing heavy atoms detected: "
        f"{n_missing_atoms}"
    )
    print(
        f"Missing terminal atoms detected: "
        f"{n_terminal_atoms}"
    )

    # Do not create artificial terminal heavy atoms at
    # crystallographic/engineered chain boundaries.
    fixer.missingTerminals = {}

    fixer.addMissingAtoms()

    fixer.addMissingHydrogens(
        pH=7.4
    )

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(OUTPUT, "w") as handle:

        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            handle,
            keepIds=True,
        )

    print()
    print("Repair complete.")
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
