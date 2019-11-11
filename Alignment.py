from config import *

""" dia(): align sequences by putting same character from both sequences """
def dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2):
    alignedSequence1 += sequence1[j - 1]
    alignedSequence2 += sequence2[i - 1]
    i -= 1
    j -= 1
    return i, j, alignedSequence1, alignedSequence2

""" left(): align sequences by adding gap in left side sequence"""
def left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2):
    alignedSequence1 += sequence1[j - 1]
    alignedSequence2 += '-'
    j -= 1
    return i, j, alignedSequence1, alignedSequence2

""" up(): align sequences by adding gap in up sequence"""
def up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2):
    alignedSequence1 += '-'
    alignedSequence2 += sequence2[i - 1]
    i -= 1
    return i, j, alignedSequence1, alignedSequence2

def align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2):
    while j > 0:
        alignedSequence1 += sequence1[j - 1]
        alignedSequence2 += '-'
        j -= 1
    while i > 0:
        alignedSequence1 += '-'
        alignedSequence2 += sequence2[i - 1]
        i -= 1
    return alignedSequence1, alignedSequence2, i, j

##########################################################################
# needlemanWunschAlign():                                                #
# 1. Traceback from bottom right cell to left top most element           #
# 2. Add gap in sequence to match highest number of symbols in sequences #
# 3. Get all the possible alignment of sequences                         #
##########################################################################

def needlemanWunschAlign(sequence1, sequence2, scoringMatrix):
    as1 = []
    as2 = []
    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            # Match current score with diagonal cell (if same then no gap is added)
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1,sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1,sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1,sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    alignedSequence1, alignedSequence2 = "", ""
    i, j = len(sequence2), len(sequence1)
    while i > 0 and j > 0:
        """ Traceback Starting from last cell of matrix to ending with top most element"""
        if scoringMatrix[i][j] == scoringMatrix[i][j - 1] + GAP_PENALTY:  # Match current score with left side cell
            i, j, alignedSequence1, alignedSequence2 = left(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        if scoringMatrix[i][j] == scoringMatrix[i - 1][j] + GAP_PENALTY:  # Match current score with up side cell
            i, j, alignedSequence1, alignedSequence2 = up(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
        elif scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF):
            """" Match current score with diagonal cell (if same then no gap is added)"""
            i, j, alignedSequence1, alignedSequence2 = dia(i, j, alignedSequence1, alignedSequence2, sequence1, sequence2)
    # Finish tracing up to the top most cell
    alignedSequence1, alignedSequence2, i, j = align(alignedSequence1, alignedSequence2, i, j, sequence1, sequence2)
    as1.append(alignedSequence1)
    as2.append(alignedSequence2)

    sequences = list(dict.fromkeys(zip(as1, as2)))
    sequences = [list(t) for t in zip(*sequences)]
    for i in range(0,len(sequences)):
        for j in range(0,len(sequences[i])):
            sequences[i][j] = sequences[i][j][::-1]

    return sequences
