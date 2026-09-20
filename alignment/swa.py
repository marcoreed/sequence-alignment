from dataclasses import dataclass
import numpy as np
from enum import IntEnum

@dataclass
class AlignmentResult:
    final_x_idx: int
    final_y_idx: int
    similarity: float
    seq_x_aligned: str
    seq_y_aligned: str

class Trace(IntEnum):
    STOP = 0
    UP = 1
    DIAGONAL= 2
    LEFT = 3

def _fill_cell(X, Y, A, row, col, match, mismatch, gap):
    similar = match if X[row - 1] == Y[col - 1] else mismatch
    diagonal_score = A[row - 1, col - 1] + similar
    up_score = A[row - 1, col] + gap
    left_score = A[row, col - 1] + gap

    best = max(0, diagonal_score, up_score, left_score)

    if best == 0:
        trace = Trace.STOP
    elif best == diagonal_score:
        trace = Trace.DIAGONAL
    elif best == up_score:
        trace = Trace.UP
    else:
        trace = Trace.LEFT

    return best, trace

def smith_waterman(X: str, Y: str, match = 1, mismatch = -1, gap = -1) -> list[list]:
    """
    Args:
        X: DNA sequence X
        Y: DNA sequence Y
        match (int): How much to reward a match
        mismatch (int): How much to punish a mismatch
        gap (int): How much to punish a gap
    """
    # Initialise 0's matrix with X vertical and Y horizontal
    row_len = len(X) + 1
    col_len = len(Y) + 1

    # Initialise the cumulative similarity matrix and trace matrix
    A = np.zeros((row_len, col_len))
    T = np.zeros((row_len, col_len))

    print('Performing forward SWA pass...')

    # For each cell (going top to bottom, left to right) carry out the SWA and record the trace
    for row in range(1, row_len):
        for col in range(1, col_len):
            A[row, col], T[row, col] = _fill_cell(X, Y, A, row, col, match, mismatch, gap)

    print('Calculating similarity...')

    # Find the index of the maximum value in A
    max_val = np.max(A)
    all_max_indices = np.argwhere(A==max_val)
    # Find the last occurrence of the highest value
    max_row, max_col = all_max_indices[-1]

    # Save the similarity score by normalising max_val by the maximum possible score
    min_seq_len = min(len(X), len(Y))
    similarity_score = max_val/min_seq_len


    # Initialise the aligned sequences
    X_aligned = ""
    Y_aligned = ""

    print('Backtracking...')

    # Backtrack the trace, saving the aligned sequence on tbe way
    while int(T[max_row, max_col]) != Trace.STOP:
        if T[max_row, max_col] == Trace.DIAGONAL:
            X_aligned = X[max_row-1] + X_aligned
            Y_aligned = Y[max_col-1] + Y_aligned
            max_row -= 1
            max_col -= 1
        elif T[max_row, max_col] == Trace.LEFT:
            X_aligned = "-" + X_aligned
            Y_aligned = Y[max_col-1] + Y_aligned
            max_col -= 1
        elif T[max_row, max_col] == Trace.UP:
            X_aligned = X[max_row-1] + X_aligned
            Y_aligned = "-" + Y_aligned
            max_row -= 1
    
    return AlignmentResult(
        similarity=similarity_score,
        final_x_idx=max_row - 1,
        final_y_idx=max_col - 1,
        seq_x_aligned=X_aligned,
        seq_y_aligned=Y_aligned
    )

if __name__ == "__main__":
    result = smith_waterman("AAABAAGDEFWAAA", "ILFFACAAAAA")
    print((result.seq_y_aligned, result.seq_x_aligned))
