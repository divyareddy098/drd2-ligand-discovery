"""Represent discontinuous 6CM4 DRD2 coordinates as separate chains.

Experimental coordinate segments:
    A: 35-139
    B: 144-222
    C: 364-442

Coordinates and residue numbering are preserved.
Only chain identifiers are changed so preparation software does not
infer peptide bonds across crystallographic gaps.
"""

from pathlib import Path


INPUT = Path(
    "structures/processed/6CM4_DRD2_receptor.pdb"
)

OUTPUT = Path(
    "structures/processed/6CM4_DRD2_split.pdb"
)


def chain_for_residue(resnum):

    if 35 <= resnum <= 139:
        return "A"

    if 144 <= resnum <= 222:
        return "B"

    if 364 <= resnum <= 442:
        return "C"

    return None


def main():

    output_lines = []

    counts = {
        "A": 0,
        "B": 0,
        "C": 0,
    }

    previous_chain = None

    with open(INPUT) as handle:

        for line in handle:

            if not line.startswith("ATOM"):
                continue

            resnum = int(
                line[22:26].strip()
            )

            new_chain = chain_for_residue(
                resnum
            )

            if new_chain is None:
                continue

            # Add TER when moving between experimental
            # coordinate segments.
            if (
                previous_chain is not None
                and new_chain != previous_chain
            ):
                output_lines.append("TER\n")

            # PDB chain ID is column 22
            new_line = (
                line[:21]
                + new_chain
                + line[22:]
            )

            output_lines.append(
                new_line
            )

            counts[new_chain] += 1
            previous_chain = new_chain

    output_lines.append("TER\n")
    output_lines.append("END\n")

    OUTPUT.write_text(
        "".join(output_lines)
    )

    print("=" * 60)
    print("6CM4 DRD2 CHAIN-BREAK PREPARATION")
    print("=" * 60)

    print(
        f"Segment A (35-139): "
        f"{counts['A']:,} atoms"
    )

    print(
        f"Segment B (144-222): "
        f"{counts['B']:,} atoms"
    )

    print(
        f"Segment C (364-442): "
        f"{counts['C']:,} atoms"
    )

    print(
        f"Total atoms: "
        f"{sum(counts.values()):,}"
    )

    print(f"\nSaved: {OUTPUT}")


if __name__ == "__main__":
    main()
