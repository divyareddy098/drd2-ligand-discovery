"""Extract human DRD2 and crystallographic risperidone from PDB 6CM4.

6CM4 contains an engineered T4 fusion protein within chain A.
According to the PDB DBREF records:

DRD2: residues 35-222
T4 fusion: residues 1002-1161
DRD2: residues 362-443

Only the human DRD2 portions are retained for receptor preparation.
"""

from pathlib import Path


INPUT = Path("structures/raw/6CM4.pdb")

OUTPUT_DIR = Path("structures/processed")

RECEPTOR = OUTPUT_DIR / "6CM4_DRD2_receptor.pdb"

LIGAND = OUTPUT_DIR / "6CM4_risperidone.pdb"

LIGAND_CODE = "8NU"

DRD2_RANGES = [
    (35, 222),
    (362, 443),
]


def is_drd2_residue(residue_number):
    """Return True if residue belongs to human DRD2."""

    return any(
        start <= residue_number <= end
        for start, end in DRD2_RANGES
    )


def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    receptor_lines = []
    ligand_lines = []

    with open(INPUT) as handle:

        for line in handle:

            if line.startswith("ATOM"):

                chain = line[21].strip()

                try:
                    residue_number = int(
                        line[22:26].strip()
                    )
                except ValueError:
                    continue

                if (
                    chain == "A"
                    and is_drd2_residue(
                        residue_number
                    )
                ):
                    receptor_lines.append(line)

            elif line.startswith("HETATM"):

                residue_name = (
                    line[17:20].strip()
                )

                if residue_name == LIGAND_CODE:
                    ligand_lines.append(line)

    receptor_lines.append("END\n")
    ligand_lines.append("END\n")

    RECEPTOR.write_text(
        "".join(receptor_lines)
    )

    LIGAND.write_text(
        "".join(ligand_lines)
    )

    print("=" * 60)
    print("6CM4 HUMAN DRD2 EXTRACTION")
    print("=" * 60)

    print(
        f"DRD2 atoms retained: "
        f"{len(receptor_lines) - 1:,}"
    )

    print(
        f"Risperidone atoms: "
        f"{len(ligand_lines) - 1:,}"
    )

    print("\nHuman DRD2 residue ranges:")
    print("35-222")
    print("362-443")

    print(
        f"\nReceptor: {RECEPTOR}"
    )

    print(
        f"Ligand:   {LIGAND}"
    )


if __name__ == "__main__":
    main()
