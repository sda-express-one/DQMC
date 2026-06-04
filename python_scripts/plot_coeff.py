#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Plot two-column data as points."
    )

    parser.add_argument(
        "input_file",
        help="Input text file with two columns: x y"
    )

    parser.add_argument(
        "-o", "--output",
        default="plot.jpg",
        help="Output image file (default: plot.jpg)"
    )

    parser.add_argument(
        "--log",
        action="store_true",
        help="Use logarithmic y-axis"
    )

    parser.add_argument(
        "--abs",
        action="store_true",
        help="Plot absolute value of y"
    )

    args = parser.parse_args()

    # Load data
    data = np.loadtxt(args.input_file)

    x = data[:, 0]
    y = data[:, 1]

    # Apply modulus if requested
    if args.abs:
        y = np.abs(y)

    # For log plots, remove non-positive values
    
    # Log plot => always use modulus
    if args.log:
        y = np.abs(y)

        # Remove zeros (cannot be shown on a log scale)
        mask = y > 0
        x = x[mask]
        y = y[mask]

        if len(y) == 0:
            raise ValueError(
                "No non-zero values remain for logarithmic plotting."
            )

    # Create figure
    plt.figure(figsize=(6, 4))

    # Points only
    plt.plot(
        x,
        y,
        linestyle="--",
        linewidth=0.5,
        marker=".",
        markersize=3,
        color="blue"
    )

    plt.xlabel("x")
    plt.ylabel("|y|" if args.abs else "y")

    if args.log:
        plt.yscale("log")

    plt.tight_layout()

    plt.savefig(args.output, dpi=300, format="jpg")
    plt.close()

    print(f"Saved plot to {args.output}")


if __name__ == "__main__":
    main()