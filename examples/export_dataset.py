"""Example 3: Building and exporting a multi-molecule benchmark dataset with provenance manifest."""

from pathlib import Path
import tempfile

from cccbdb import CCCBDB

FIXTURES = Path(__file__).resolve().parent.parent / "tests" / "fixtures"


def main() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir)
        with CCCBDB(cache=True, cache_dir=out_dir / "cache", offline=True, use_sqlite=True) as db:
            # Seed cache with experimental fixtures for H2O and CO
            db.disk_cache.set(
                db.transport.build_url("exp2x.asp", {"casno": "7732185"}),
                (FIXTURES / "exp2x_h2o.html").read_text(encoding="utf-8"),
            )
            db.disk_cache.set(
                db.transport.build_url("exp2x.asp", {"casno": "630080"}),
                (FIXTURES / "exp2x_co.html").read_text(encoding="utf-8"),
            )

            dataset = db.collect(
                molecules=["H2O", "CO"],
                properties=["dipole", "zpe", "enthalpy"],
                method="Experimental",
                basis_set=None,
                delay=0.0,
            )

            csv_file = out_dir / "experimental_benchmark.csv"
            dataset.to_csv(csv_file)
            print("=== Exported CSV Content ===")
            print(csv_file.read_text(encoding="utf-8"))

            manifest_file = Path(str(csv_file) + ".manifest.json")
            print("=== Exported Reproducibility Manifest ===")
            print(manifest_file.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
