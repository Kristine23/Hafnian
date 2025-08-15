#!/usr/bin/env python3.10

import argparse
import numpy as np
from hafnian import HafnianCalculator


def main():
    parser = argparse.ArgumentParser(
    description="Compute the hafnian of a matrix from a file, "
                    "or use a 6x6 block-diagonal matrix with blocks [[0,2],[2,0]] "
                    "(which has exactly one perfect matching)."
    )
    parser.add_argument(
        "matrix_file", nargs="?", default=None, 
        help="Path to .npy or .csv matrix file (optional)"
    )
    parser.add_argument(
        "--blocks", type=int, default=1,
        help="Number of blocks [[0,2][2,0]] on default block-diagonal matrix (must be >=1). Used if no matrix_file is provided. Default: 2"
    )
    args = parser.parse_args()


    if args.matrix_file:
        # Load from file
        try:
            if args.matrix_file.endswith('.npy'):
                matrix = np.load(args.matrix_file)
            else:
                matrix = np.loadtxt(args.matrix_file, delimiter=',')
        except Exception as e:
            print("Failed to load file:", e)
            return
    else:
        # Use a block-diagonal matrix
        matrix = np.kron(np.identity(args.blocks), [[0, 2], [2, 0]])

    result = HafnianCalculator.hafnian(matrix)
    print("Hafnian:", result)
    print("VERY rough estimate of the numerical precision", HafnianCalculator.estimate_precision(matrix))

if __name__ == '__main__':
    main()