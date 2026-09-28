"""Calculate docking-box center from crystallographic risperidone."""

from pathlib import Path

import numpy as np


LIGAND = Path(
    "structures/processed/6CM4_risperidone.pdb"
)


def main():

    coordinates = []

    with open(LIGAND) as handle:

        for line in handle:

            if line.startswith(
                ("ATOM", "HETATM")
            ):

                x = float(
                    line[30:38]
                )

                y = float(
                    line[38:46]
                )

                z = float(
                    line[46:54]
                )

                coordinates.append(
                    [x, y, z]
                )

    coordinates = np.asarray(
        coordinates
    )

    center = coordinates.mean(
        axis=0
    )

    minimum = coordinates.min(
        axis=0
    )

    maximum = coordinates.max(
        axis=0
    )

    dimensions = (
        maximum - minimum
    )

    print("=" * 55)
    print("6CM4 BINDING SITE")
    print("=" * 55)

    print(
        f"Ligand atoms: "
        f"{len(coordinates)}"
    )

    print("\nLigand centroid")

    print(
        f"center_x = {center[0]:.3f}"
    )

    print(
        f"center_y = {center[1]:.3f}"
    )

    print(
        f"center_z = {center[2]:.3f}"
    )

    print("\nLigand dimensions")

    print(
        f"x = {dimensions[0]:.3f} Å"
    )

    print(
        f"y = {dimensions[1]:.3f} Å"
    )

    print(
        f"z = {dimensions[2]:.3f} Å"
    )


if __name__ == "__main__":
    main()
