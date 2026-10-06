"""Check release-part hashes and restore the exact GTEx v10 validation inputs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from .paper import ROOT, sha256


def restore(bundle_dir: Path, destination: Path) -> None:
    manifest = json.loads((ROOT / "reproducibility/input_bundle.json").read_text())
    destination.mkdir(parents=True, exist_ok=True)
    for filename, entry in manifest["files"].items():
        target = destination / filename
        if target.exists():
            if sha256(target) == entry["sha256"]:
                print(filename + ": already restored and verified.")
                continue
            raise FileExistsError(f"Refusing to overwrite a different existing input: {target}")
        temporary = destination / (filename + ".restoring")
        if temporary.exists():
            raise FileExistsError(f"A previous restore left {temporary}; inspect it before retrying.")
        try:
            with temporary.open("xb") as output:
                for part in entry["parts"]:
                    source = bundle_dir / part["name"]
                    if not source.is_file() or source.stat().st_size != part["bytes"] or sha256(source) != part["sha256"]:
                        raise RuntimeError(f"Missing or modified bundle part: {source}")
                    with source.open("rb") as stream:
                        shutil.copyfileobj(stream, output, length=8 * 1024 * 1024)
            if temporary.stat().st_size != entry["bytes"] or sha256(temporary) != entry["sha256"]:
                raise RuntimeError("Restored input checksum mismatch: " + filename)
            temporary.replace(target)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
        print(filename + ": restored and SHA-256 verified.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--destination", type=Path, default=ROOT / "data/parquet/v10")
    args = parser.parse_args()
    restore(args.bundle_dir.resolve(), args.destination.resolve())


if __name__ == "__main__":
    main()
