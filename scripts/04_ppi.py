"""Integrate STRING PPI network with switching gene sets."""
import argparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", choices=["v8", "v10"], required=True)
    args = parser.parse_args()
    print(f"TODO: PPI integration for GTEx {args.version}")


if __name__ == "__main__":
    main()
