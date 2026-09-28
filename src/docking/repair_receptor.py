"""Repair missing atoms in the experimental 6CM4 DRD2 receptor."""

from pathlib import Path

from pdbfixer import PDBFixer
from openmm.app import PDBFile


INPUT = Path(
    "structures/processed/6CM4_DRD2_receptor.pdb"
)

OUTPUT = Path(
    "structures/processed/6CM4_DRD2_receptor_repaired.pdb"
)


def main():

    print("=" * 60)
    print("6CM4 DRD2 RECEPTOR REPAIR")
    print("=" * 60)
    print(f"Input: {INPUT}")

    fixer = PDBFixer(
        filename=str(INPUT)
    )

    # Identify nonstandard residues.
    fixer.findNonstandardResidues()

    print(
        f"Nonstandard residues: "
        f"{len(fixer.nonstandardResidues)}"
    )

    # PDBFixer requires missing residues to be initialized
    # before findMissingAtoms().
    fixer.findMissingResidues()

    detected_missing_residues = dict(
        fixer.missingResidues
    )

    print(
        f"Missing residue segments detected: "
        f"{len(detected_missing_residues)}"
    )

    # Do not reconstruct completely unresolved residues.
    fixer.missingResidues = {}

    # Detect missing atoms in residues that are present.
    fixer.findMissingAtoms()

    n_missing_atoms = sum(
        len(atoms)
        for atoms in fixer.missingAtoms.values()
    )

    n_missing_terminals = sum(
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
        f"{n_missing_terminals}"
    )

    # IMPORTANT:
    # Do not construct terminal heavy atoms at artificial
    # crystallographic chain boundaries.
    fixer.missingTerminals = {}

    # Repair missing side-chain/heavy atoms.
    fixer.addMissingAtoms()

    # Add hydrogens at approximately physiological pH.
    fixer.addMissingHydrogens(
        pH=7.4
    )

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
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
