"""Run the full STAMP pipeline (steps 01 through 06) for one GTEx version."""
import argparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", choices=["v8", "v10"], required=True)
    args = parser.parse_args()
    print(f"TODO: run full pipeline for GTEx {args.version}")


if __name__ == "__main__":
    main()
