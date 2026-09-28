"""Download the experimental DRD2 structure from RCSB PDB."""

from pathlib import Path

import requests


PDB_ID = "6CM4"

URL = (
    f"https://files.rcsb.org/download/"
    f"{PDB_ID}.pdb"
)

OUTPUT = Path(
    f"structures/raw/{PDB_ID}.pdb"
)


def main():

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    print(
        f"Downloading PDB {PDB_ID}..."
    )

    response = requests.get(
        URL,
        timeout=60
    )

    response.raise_for_status()

    OUTPUT.write_bytes(
        response.content
    )

    print(
        f"Saved: {OUTPUT}"
    )

    print(
        f"File size: "
        f"{OUTPUT.stat().st_size / 1024:.1f} KB"
    )


if __name__ == "__main__":
    main()
