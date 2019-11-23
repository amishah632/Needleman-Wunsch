from config import *
import re
from sys import *
from numpy import *
import argparse as ap

opfile = ""
cnt = 1
aseq1, match_sym, aseq2 = [], [], []

def traceback(i, j, left_path, up_path, diagonal_path, sequence1, sequence2, aligned_seq1, aligned_seq2):
    """
    traceback whole matrix from bottom-right element to top-left element
    to get all possible alignments
    :param i: row index
    :param j: column index
    :param left_path: matrix that directs to left path
    :param up_path: matrix that directs to up path
    :param diagonal_path: matrix that directs to diagonal path
    :param sequence1: sequence 1
    :param sequence2: sequence 2
    :param aligned_seq1: calculated alignment for sequence 1
    :param aligned_seq2: calculated alignment for sequence 2
    """
    match_symbol = []
    global cnt
    if cnt <= MAX_NUMBER_PATHS:
        if i == j == 0:
            aligned_seq1.append(sequence1[i])
            aligned_seq2.append(sequence2[j])
            aligned_seq1, aligned_seq2 = aligned_seq1[::-1], aligned_seq2[::-1]
            for symbol in range(len(aligned_seq1)):
                match_symbol.append('|') if aligned_seq1[symbol] == aligned_seq2[symbol] else match_symbol.append(' ')
            as1, match, as2 = ''.join(aligned_seq1), ''.join(match_symbol), ''.join(aligned_seq2)
            cnt += 1
            global aseq1, aseq2, match_sym
            aseq1.append(as1)
            aseq2.append(as2)
            match_sym.append(match)
        else:
            if diagonal_path[i, j] == 1:
                temp_seq1,temp_seq2 = aligned_seq1.copy(), aligned_seq2.copy()
                temp_seq1.append(sequence1[i])
                temp_seq2.append(sequence2[j])
                traceback(i - 1, j - 1, left_path, up_path, diagonal_path, sequence1, sequence2, temp_seq1, temp_seq2)
            if left_path[i, j] == 1:
                temp_seq1,temp_seq2 = aligned_seq1.copy(), aligned_seq2.copy()
                temp_seq1.append('-')
                temp_seq2.append(sequence2[j])
                traceback(i, j - 1, left_path, up_path, diagonal_path, sequence1, sequence2, temp_seq1, temp_seq2)
            if up_path[i, j] == 1:
                temp_seq1,temp_seq2 = aligned_seq1.copy(), aligned_seq2.copy()
                temp_seq1.append(sequence1[i])
                temp_seq2.append('-')
                traceback(i - 1, j, left_path, up_path, diagonal_path, sequence1, sequence2, temp_seq1, temp_seq2)

def align_sequences(sequence1, sequence2):
    """
    :param sequence1: sequence 1
    :param sequence2: sequence 2
    :return: scoring_matrix, alignments for both sequences and indicates matching elements from both sequences
    """

    scoring_matrix = zeros((len(sequence1) + 1, len(sequence2) + 1), int)
    left_path, up_path, diagonal_path = (zeros(shape=(len(sequence1), len(sequence2))) for _ in range(3))

    """ compute matrix elements and fill into matrix"""
    scoring_matrix[0][0] = 0  # Set 0th row 0th column element to zero
    for i in range(1, len(sequence1) + 1):
        scoring_matrix[i][0] = scoring_matrix[i - 1][0] + GAP_PENALTY  # Fill 0th column by adding Gap penalty to each cell
    for j in range(1, len(sequence2) + 1):
        scoring_matrix[0][j] = scoring_matrix[0][j - 1] + GAP_PENALTY  # Fill 0th row by adding Gap penalty to each cell
    for i in range(1, len(sequence1) + 1):
        for j in range(1, len(sequence2) + 1):
            value_derived_from_up = scoring_matrix[i - 1][j] + GAP_PENALTY
            value_derived_from_left = scoring_matrix[i][j - 1] + GAP_PENALTY
            value_derived_from_diagonal = \
                scoring_matrix[i - 1][j - 1] + (SAME if sequence1[i - 1] == sequence2[j - 1] else DIFF)
            scoring_matrix[i][j] = max(value_derived_from_up, value_derived_from_left, value_derived_from_diagonal)
            if scoring_matrix[i][j] == value_derived_from_left:
                left_path[i - 1, j - 1] = 1  # denote left path matri
            if scoring_matrix[i][j] == value_derived_from_up:
                up_path[i - 1, j - 1] = 1  # denote up path matrix
            if scoring_matrix[i][j] == value_derived_from_diagonal:
                diagonal_path[i - 1, j - 1] = 1  # denote diagonal path matrix
    traceback(len(sequence1) - 1, len(sequence2) - 1, left_path, up_path, diagonal_path, sequence1, sequence2, [], [])
    global aseq1, aseq2, match_sym
    return scoring_matrix, aseq1, aseq2, match_sym

if __name__ == "__main__":
    """ Values taken from Command Line Arguments """
    parser = ap.ArgumentParser()
    parser.add_argument("sequence1", help="Sequence1 FASTA File path")
    parser.add_argument("sequence2", help="Sequence2 FASTA File path")
    parser.add_argument("outputFile", help="Output FileName")
    parser.add_argument("--gapPenalty", help="Gap Penalty Point (-ve less than 0)")
    parser.add_argument("--match", help="Match Point (+ve greater than 0)")
    parser.add_argument("--mismatch", help="Mismatch Point (-ve less than 0)")
    args = parser.parse_args()

    if args.gapPenalty:
        """" If user has given gap penalty in command line argument then replace it as GAP_PENALTY """
        if int(args.gapPenalty) >= 0:
            print("Gap penalty must be negative integer")
            exit(1)
        else:
            GAP_PENALTY = int(args.gapPenalty)
    if args.match:
        """" If user has given match points in command line argument then replace it as SAME """
        if int(args.match) <= 0:
            print("Match Point must be positive integer")
            exit(1)
        else:
            SAME = int(args.match)
    if args.mismatch:
        """" If user has given mismatch points in command line argument then replace it as DIFF """
        if int(args.mismatch) >= 0:
            print("Mismatch point must be negative integer")
            exit(1)
        else:
            DIFF = int(args.mismatch)

    s1_fasta = open(args.sequence1, "r")
    s2_fasta = open(args.sequence2, "r")
    opfile = open(args.outputFile, "w+")

    sequence1,sequence2 = s1_fasta.read(),s2_fasta.read()

    if sequence1[0] == '>':
        sequence1 = sequence1.partition("\n")
        sequence1 = sequence1[2].replace("\n", "").replace(" ", "").upper()
    else:
        sequence1 = sequence1.replace("\n", "").replace(" ", "").upper()

    if sequence2[0] == '>':
        sequence2 = sequence2.partition("\n")
        sequence2 = sequence2[2].replace("\n", "").replace(" ", "").upper()
    else:
        sequence2 = sequence2.replace("\n", "").replace(" ", "").upper()

    if (len(sequence1) > MAX_SEQ_LENGTH) and (len(sequence2) > MAX_SEQ_LENGTH):
        stderr.write("Length of both sequences are too long")
        exit(1)
    elif len(sequence1) > MAX_SEQ_LENGTH:
        stderr.write("Length of sequence1 is too long")
        exit(1)
    elif len(sequence2) > MAX_SEQ_LENGTH:
        stderr.write("Length of sequence2 is too long")
        exit(1)
    else:
        if bool(re.match('^[GCAUT]+$', sequence1)) and bool(re.match('^[GCAUT]+$', sequence2)):
            scoring_matrix, as1, as2, match = align_sequences(sequence1, sequence2)
            print("\nScoring Matrix: \n", scoring_matrix, "\n\nScore: ", scoring_matrix[len(sequence1)][len(sequence2)])
            for i in range(0, len(as1)):
                opfile.write("Alignment: %s\n%s\n%s\n%s\n\n" % (i+1, as1[i], match[i], as2[i]))
                print("\nAlignment ", i+1, ":\n", as1[i], "\n", match[i], "\n", as2[i])
            opfile.write("Score: %s" % (scoring_matrix[len(sequence1)][len(sequence2)]))
            opfile.close()
        else:
            print("Only DNA or mRNA sequence is allowed")
            exit(1)
