#!/usr/bin/env python3.10
import numpy as np

class HafnianCalculator:
    
    @staticmethod
    def _powerset(int_set: list[int]) -> list[list[int]]:
        """
        Generate all possible subsets (the power set) of the given list of integers .
        """
        
        n = len(int_set)
        powerset = []
        # For each binary number in range 0 to 2^n - 1, use its bits to select elements
        for i in range(1 << n):
            powerset.append([int_set[j] for j in range(n) if (i & (1 << j))])
        return powerset
    
    @staticmethod
    def _column_swapped(matrix: np.ndarray) -> np.ndarray:
        """
        Swap adjacent columns in the given square matrix in 2-by-2 blocks.
        
        For a matrix of size n x n, where n is even, columns 0 and 1 are swapped,
        columns 2 and 3 are swapped, etc.
        """
        
        n = len(matrix)
        half_n = n // 2
        # Construct a block-diagonal matrix that swaps each pair of columns
        column_swapping_matrix  = np.kron(np.identity(half_n), [[0, 1], [1, 0]])
        # Perform the column swaps by matrix multiplication
        swapped_matrix= np.matmul(column_swapping_matrix, matrix)
        return swapped_matrix
    
    @staticmethod
    def _principal_submatrix(matrix: np.ndarray, indices: list[int]) -> np.ndarray:
        """
        Extract the principal submatrix of a square matrix corresponding to a set of indices.
        """
        submatrix = matrix[np.ix_(indices, indices)]
        return submatrix
    
    @staticmethod
    def _list_inner_sum_coefficients(matrix: np.ndarray, half_n: int) -> list[float]:
        """
        For each integer k in 1..n//2, compute the sum of eigenvalues raised to the k-th power
        and divide by 2k, i.e.:

            coefficient_k = sum_i (lambda_i ** k) / (2 * k)

        where lambda_i are the eigenvalues of the input matrix.

        Returns:
            list of float: The computed coefficients for each k (k=1 to n//2).
        """
            
        eigenvalues = np.linalg.eigvals(matrix)
        dim = len(eigenvalues)
        power_trace = np.ones(dim)
        
        inner_coefficient_list = []
        for k in range(1, half_n+1):
            
            power_trace = np.array([power_trace[i]*eigenvalues[i] for i in range(dim)])  # eigenvalues ** k
            inner_coefficient = sum(power_trace / (2*k))
            inner_coefficient_list.append(inner_coefficient)
        
        return inner_coefficient_list
        
    @staticmethod
    def _coefficient_in_generating_function(matrix: np.ndarray, half_n: int) -> float:
        """
        Compute the coefficient of order n/2 (=half_n) in the expansion of the generating function (Equation 3.21):

            p(λ, B) = Σ_{j=1}^{n/2} (1/j!) * (Σ_{k=1}^{n/2} (tr(B^k) * λ^k)/(2k))^j

        This function computes the coefficient of λ^{half_n} in the expansion of p(λ, B). 

        In this implementation, the terms
            (tr(B^k) * λ^k) / (2k)
        are denoted as the inner_coefficients.

        """

        inner_coefficients = HafnianCalculator._list_inner_sum_coefficients(matrix, half_n)
        scalar = 1.0

        # Stores coefficients of the generating function for λ^i (i = 1...half_n) for the previous j-1 value
        previous_inner_coefficients = inner_coefficients[:]

        # Initialize result with the λ^{half_n} coefficient for j=1
        result_coefficient = inner_coefficients[half_n - 1]

        # Iteratively build higher-order coefficients (for j=2 to half_n) using convolution-like accumulation
        for j in range(2, half_n + 1):
            current_inner_coefficients = np.zeros(half_n)
            for d1 in range(1, half_n + 1):
                for d2 in range(1, half_n + 1): 
                    idx = d1 + d2 - 1
                    if idx < half_n:
                        # Combine products of pairs of coefficients to build higher order
                        current_inner_coefficients[idx] += (
                            inner_coefficients[d2-1] * previous_inner_coefficients[d1-1]
                        )
            previous_inner_coefficients = current_inner_coefficients[:]
            scalar /= j
            # Aggregate the coefficient of λ^{half_n} for the current j
            result_coefficient += scalar * current_inner_coefficients[half_n - 1]

        return result_coefficient
        
    @staticmethod    
    def _check_input(matrix_A: np.ndarray) -> None:
        """
        Check that the input matrix is square, has even dimensions, and is symmetric.
        """
        shape = matrix_A.shape

        # Check matrix is square
        if shape[0] != shape[1]:
            raise ValueError("The matrix is not square.")
        
        # Check matrix dimension is even
        if shape[0] % 2 != 0:
            raise ValueError("The size of the matrix is not even.")
        
        # Check matrix is symmetric
        if not np.allclose(matrix_A, matrix_A.T, rtol=1e-12, atol=1e-12):
            raise ValueError("The matrix is not symmetric.")
        
        # Check the matrix is 2D
        if matrix_A.ndim != 2:
            raise ValueError("The matrix is not 2D")
        
    @staticmethod
    def hafnian(matrix_A: np.ndarray) -> float:
        """
        Compute the hafnian of a symmetric, even-dimensional matrix.

        Uses Equation (3.24) to calculate the hafnian:

            haf(A) = Σ_{Z ∈ P([n/2])} (−1)^{n/2−|Z|} * (coefficient of λ^{n/2} in p(λ, B))

        where:
            - P([n/2]) is the powerset of {0, ..., n/2 - 1},
            - |Z| is the cardinality of subset Z,
            - p(λ, B) is the generating function (Equation 3.21):

                p(λ, B) = Σ_{j=1}^{n/2} (1/j!) * (Σ_{k=1}^{n/2} (tr(B^k) * λ^k) / (2k) )^j
        """

        HafnianCalculator._check_input(matrix_A)

        n = len(matrix_A)
        half_n = n // 2

        # Create all subsets of {0, 1, ..., half_n - 1}
        integer_set = list(range(half_n))
        powerset = HafnianCalculator._powerset(integer_set)

        # Swap columns in pairs 
        column_swapped_A = HafnianCalculator._column_swapped(matrix_A)

        hafnian = 0.0
        # Sum over all subsets in the powerset
        for subset in powerset:
            set_cardinality = len(subset)
            sign = (-1) ** (half_n - set_cardinality)

            # Map the subset to row/column indices (expanded for each pair)
            row_col_idx = [2*i + o for i in subset for o in (0, 1)]
            principal_submatrix = HafnianCalculator._principal_submatrix(column_swapped_A, row_col_idx)

            # Compute the coefficient for λ^{half_n}
            coefficient = HafnianCalculator._coefficient_in_generating_function(principal_submatrix, half_n)

            hafnian += sign * coefficient
        
        return hafnian
        
    @staticmethod
    def estimate_precision(matrix: np.ndarray) -> float:
        
        """
        VERY rough estimate of the numerical precision for the implemented eigenvalue-based calculations
        on a square matrix. This estimate accounts for error amplification arising from repeated
        matrix multiplications and combinatorial sums (as encountered in the hafnian calculations).
        """
        
        n = len(matrix)
        eigenvalues = np.linalg.eigvals(matrix)
        
        # As per floats in Python and LAPACK/NumPy documentation for eigenvalues. 
        base_precision = 1e-15
        
        # Use the absolute value of the largest eigenvalue for safety 
        # (works for both real symmetric and complex Hermitian matrices).
        # By the interlacing theorem, every principal submatrix has 
        # its largest eigenvalue ≤ the largest eigenvalue of the original Hermitian matrix.
        max_eigenvalue = np.max(np.abs(eigenvalues))
        
        # Estimate error amplification: Power of eigenvalue 
        amplification = max_eigenvalue ** (n // 2)
        
        # Estimate based on the sum of powers of the eigenvalues
        estimated = amplification * base_precision * (n / 2)
        
        #Multipyer due to summing over all powersets
        total_estimated_precision = estimated * (2 ** n)
        
        return total_estimated_precision


        



