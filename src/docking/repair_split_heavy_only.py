"""Repair missing heavy atoms in split 6CM4 DRD2 without adding hydrogens."""

from pathlib import Path

from pdbfixer import PDBFixer
from openmm.app import PDBFile


INPUT = Path(
    "structures/processed/6CM4_DRD2_split.pdb"
)

OUTPUT = Path(
    "structures/processed/6CM4_DRD2_heavy_repaired.pdb"
)


def main():

    print("=" * 60)
    print("6CM4 DRD2 HEAVY-ATOM REPAIR")
    print("=" * 60)

    fixer = PDBFixer(filename=str(INPUT))

    fixer.findNonstandardResidues()
    fixer.findMissingResidues()

    # Do not reconstruct completely missing residues.
    fixer.missingResidues = {}

    fixer.findMissingAtoms()

    n_missing = sum(
        len(atoms)
        for atoms in fixer.missingAtoms.values()
    )

    n_terminal = sum(
        len(atoms)
        for atoms in fixer.missingTerminals.values()
    )

    print(f"Missing heavy atoms: {n_missing}")
    print(f"Terminal atoms detected: {n_terminal}")

    # Do not fill artificial chain termini.
    fixer.missingTerminals = {}

    fixer.addMissingAtoms()

    # IMPORTANT:
    # Do NOT call addMissingHydrogens() here.

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

    print(f"\nSaved: {OUTPUT}")


if __name__ == "__main__":
    main()
