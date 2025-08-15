Hafnian Calculator
==================

A Python implementation to compute the hafnian of a real symmetric matrix of even size.
The implemented algorithm for computing the hafnian of the matrix is described in the following paper:
https://arxiv.org/abs/1805.12498


Comments
-----------------

- Memory usage is not efficient: the powerset is calculated and stored explicitly, requiring O(2^n) memory for an n x n matrix.
- No parallelization: the current implementation does not use parallel processing.
- Precision estimate is very rough: the provided precision estimate is only a guideline and may be wrong.
- Todo: Consider comparing to already existing implementations.


*IMPORTANT*
> This project **requires Python 3.9 or higher**.
> Older versions are not compatible due to type hinting and newer Python features.
> You can check your Python version with:
 >  python --version



Project Structure
-----------------
```
your_project/
├── README.txt
├── src/
│   ├── hafnian.py                 (HafnianCalculator class and logic)
│   ├── main.py                    (Command-line interface)
│   ├── block_diagonal.csv         (Example block-diagonal input matrix, CSV format)
│   ├── block_diagonal.npy         (Example block-diagonal input matrix, NumPy format)
└── tests/
    └── test_hafnian.py            (Unit tests)
```

Features
--------

- Compute hafnian for .npy (NumPy) or .csv matrix files.
- Uses a default block-diagonal matrix if no input file is specified.
- Includes example input files (`block_diagonal.csv` and `block_diagonal.npy`) in src/.
- Provides unittests for verification.


Getting Started
---------------

1. Make sure you are using Python 3.9 or higher.

   Check your Python version with:

       python --version

2. Install dependencies

   Use the requirements file to install all necessary packages:

       pip install -r requirements.txt

3. Run from source

   Change to the src/ directory and run:

       python main.py block_diagonal.npy

   or for CSV:

       python main.py block_diagonal.csv

   You can also use your own matrix files in .csv or .npy format.

   If no file is specified, a default block-diagonal matrix [[0,2],[2,0]] is used:

       python main.py

   You can also customize the size of the generated default matrix with the `--blocks` argument.
   Each block is [[0,2],[2,0]], so the resulting matrix will have size (2 × blocks):

       python main.py --blocks 4

   (This example generates an 8x8 block-diagonal matrix with 4 blocks.)

   The --blocks argument is ignored if you provide a matrix file.


File Formats
------------
Matrix must be real, symmetric, and of even size.

Included Example Files
----------------------

- block_diagonal.csv — Example block-diagonal matrix in CSV.
- block_diagonal.npy — Example block-diagonal matrix in NumPy format.

Run the Unit Tests
------------------

From the project root, run:

   python -m unittest discover tests

Or, if necessary:

   PYTHONPATH=src python -m unittest discover tests

Example
-------

Sample file block_diagonal.csv (found in src/):

0,1,0,0
1,0,0,0
0,0,0,1
0,0,1,0

Command to compute:

   python src/compute_hafnian.py src/block_diagonal.csv

Author
------

Kristine Vitting Klinkby Knudsen
