# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    # Creating an instance of the class
    nw = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', float(-10), float(-1))
    
    # The actual align matrix
    true_align_mat = np.array([[  0., -11., -12., -13., -14.],
                            [-11.,   5.,  -6.,  -7.,  -8.],
                            [-12.,   4.,   4.,  -1.,  -6.],
                            [-13.,  -7.,   3.,   5.,   4.]])
    
    # The actual traceback matrix
    data = [
    [None, (0, 0), (0, 1), (0, 2), (0, 3)],
    [(0, 0), (0, 0), (1, 1), (1, 2), (1, 3)],
    [(1, 0), (1, 1), (1, 1), (1, 2), (1, 3)],
    [(2, 0), (2, 1), (2, 2), (2, 2), (2, 3)]
            ]

    true_backtrace_mat = np.array(data, dtype='object')


    # Calculate the matrices from my code implementation
    my_align_mat, my_backtrace = nw.return_matrices(seq1, seq2)


    # Assert that align matrix matches the real matrix
    assert np.array_equal(true_align_mat, my_align_mat)


    # Assert that the backtrace matrix matches the real backtrace
    assert np.array_equal(true_backtrace_mat, my_backtrace)

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    

    # Creating an instance of the nw class
    nw = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', float(-10), float(-1))


    # Getting outputs of the NW align
    output = nw.align(seq3, seq4)
    score, seqA, seqB = output
    
    # Asserting that the score matches the expected score
    assert score == 17

    # Asserting the output seqA/SeqB matches the expected output SeqA/SeqB respectively
    assert seqA == 'MAVHQLIRRP'
    assert seqB == 'M---QLIRHP'




