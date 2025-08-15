import unittest
from src.hafnian import HafnianCalculator
import numpy as np
import math


class TestHafnianCalculator(unittest.TestCase):
    
    def test_unique_perfect_matching(self):
        A = np.kron(np.identity(4), [[0, 1], [1, 0]])
        result = HafnianCalculator.hafnian(A)
        #one perfect_matching
        calc_result = 1. 
        self.assertAlmostEqual(result, calc_result)
    
    def test_complete_graph(self):
        size = 6
        A = np.ones((size,size))
        np.fill_diagonal(A, 0)
        factorial = math.factorial(size)
        result = HafnianCalculator.hafnian(A)
        #number perfect matchings in a complete graph
        calc_result = math.factorial(size) / (2**(size//2) * math.factorial(size//2)) 
        self.assertAlmostEqual(result, calc_result)
    
    def test_block_matrix(self):
        number_blocks = 4
        A = np.kron(np.identity(number_blocks), [[0, 2], [2, 0]])
        result = HafnianCalculator.hafnian(A)
        #hafnian of a block-diagonal matrix is the product of the hafnians of its blocks
        calc_result = 2**number_blocks
        self.assertAlmostEqual(result, calc_result)
        
    
    def test_coefficient_in_generating_function(self):
        size = 4
        A = np.identity(size)
        result = HafnianCalculator._coefficient_in_generating_function(A, size // 2 )
        cal_result = 3
        self.assertAlmostEqual(result, cal_result)
        
    def test_non_symmetric_raises(self):
        # Should raise ValueError for non-symmetric matrices
        A = np.array([[1,2],[3,4]])
        with self.assertRaises(ValueError):
            HafnianCalculator.hafnian(A)

if __name__ == '__main__':
    unittest.main()
    




