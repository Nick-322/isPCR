#!/usr/bin/env python3

import numpy as np
def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    n = len(seq_a) + 1
    m = len(seq_b) + 1
    matrix = np.zeros((n, m))
    #####STEP 1: Fill in the first row & column for gap penalty (Initializatio)#####
    matrix[:,0] = np.linspace(0, n * gap, n)
    matrix[0,:] = np.linspace(0, m * gap, m)
    #####STEP 2: # Fill in the score matrix (Matrix-Filling)#####
    for i in range(1, n):
        for j in range(1, m):
            if seq_a[i-1] == seq_b[j-1]:
                diagonal_score = matrix[i-1][j-1] + match
            else:
                diagonal_score = matrix[i-1][j-1] + mismatch
                #choose the highest score among diagonal, vertical, and horizontal scoring
            matrix[i][j] = max(diagonal_score, matrix[i][j-1]+gap, matrix[i-1][j]+gap)
    ######STEP 3: Back-Tracing#####
    aligned_seq1 = []
    aligned_seq2 = []
    x = n - 1
    y = m - 1
    while x > 0 or y > 0:
        #vertical move = gap in horizontal sequence
        if x > 0 and matrix[x - 1][y] + gap == matrix[x][y]:
            aligned_seq1.append(seq_a[x - 1])
            aligned_seq2.append("-")
            x -= 1
        #horizontal move = gap in vertical sequence
        elif y > 0 and matrix[x][y - 1] + gap == matrix[x][y]:
            aligned_seq1.append("-")
            aligned_seq2.append(seq_b[y - 1])
            y -= 1
        #diagonal move
        else:
            aligned_seq1.append(seq_a[x - 1])
            aligned_seq2.append(seq_b[y - 1])
            x -= 1
            y -= 1
    while x > 0:
        aligned_seq1.append(seq_a[x - 1])
        aligned_seq2.append("-")
        x -= 1

    while y > 0:
        aligned_seq1.append("-")
        aligned_seq2.append(seq_b[y - 1])
        y -= 1
        
    alignment1 = ''.join(aligned_seq1)[::-1]
    alignment2 = ''.join(aligned_seq2)[::-1]
    #STEP 4: Alignment scoring
    score = 0
    for i in range(len(alignment1)):
        char1 = alignment1[i]
        char2 = alignment2[i]
        if char1 == '-' or char2 == '-':
            score += gap
        elif char1 == char2:
            score += match
        else:
            score += mismatch
    return (alignment1, alignment2), score
