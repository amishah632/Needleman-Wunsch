from config import *
import Alignment
from sys import *
from numpy import *
import argparse as ap

###############################################################################
# finalize():                                                                 #
#   Match symbols from both sequences, compute score and gives best alignment #
###############################################################################

def finalize(sequences):
    # # sequences contains tuples of form [(possibleAlignmentSeqs1),(possibleAlignmentSeqs2)]

    alignedSeq1 = sequences[0]  # sequences[0] has all possible alignments of sequence1
    alignedSeq2 = sequences[1]  # sequences[1] has all possible alignments of sequence2
    sym=[]
    scoreArray = []

    """
    Fetch every sequence one by one and match symbols of both sequences with '|' 
    if match found put '|' otherwise put ' ' 
    if match found SAME points added, if not found then DIFF points added to score 
    if gap is there then add GAP_PENALTY to score 
    """

    for x in range(0, len(alignedSeq1)):
        alignedSequence1 = alignedSeq1[x]
        alignedSequence2 = alignedSeq2[x]
        score = 0
        symbol = ""

        for i in range(0, len(alignedSequence1)):
            if alignedSequence1[i] == alignedSequence2[i]:
                symbol = symbol + '|'
                score += (SAME if alignedSequence1[i] == alignedSequence2[i] else DIFF)
            elif alignedSequence1[i] != alignedSequence2[i] and alignedSequence1[i] != '-' and alignedSequence2[i] != '-':
                score += (SAME if alignedSequence1[i] == alignedSequence2[i] else DIFF)
                symbol += ' '
            elif alignedSequence1[i] == '-' or alignedSequence2[i] == '-':
                symbol += ' '
                score += GAP_PENALTY
        sym.append(symbol)
        scoreArray.append(score)

    while len(alignedSeq1) > MAX_NUMBER_PATHS:
        alignedSeq1.pop()
        alignedSeq2.pop()
        scoreArray.pop()
        sym.pop()
    # for i in range(0, len(alignedSeq1)):
    #     print("\n\n")
    #     print(alignedSeq1[i])
    #     print(sym[i])
    #     print(alignedSeq2[i])
    #     print("Score: ",scoreArray[i])

    return sym, scoreArray, alignedSeq1, alignedSeq2

#####################################################################
# buildMatrix(): Build Matrix, Compute Matrix Elements, Fill Matrix #
#####################################################################
def buildMatrix(sequence1, sequence2):
    # rows, cols variable to set size of matrix according to length of sequence
    rows, cols = len(sequence2) + 1, len(sequence1) + 1
    scoringMatrix = zeros((rows,cols),int)

    """ compute matrix elements and fill into matrix"""
    scoringMatrix[0][0] = 0  # Set 0th row 0th column element to zero
    for i in range(1,rows):
        # Fill 0th column by adding Gap penalty to each cell
        scoringMatrix[i][0] = scoringMatrix[i-1][0] + GAP_PENALTY
    for j in range(1,cols):
        # Fill 0th row by adding Gap penalty to each cell
        scoringMatrix[0][j] = scoringMatrix[0][j-1] + GAP_PENALTY
    for i in range(1, rows):
        for j in range(1, cols):
            matrixValues = [scoringMatrix[i - 1][j] + GAP_PENALTY,  # Value derived from Up cell
                            scoringMatrix[i][j - 1] + GAP_PENALTY,  # Value derived from Left cell
                            scoringMatrix[i-1][j-1] + (SAME if sequence1[j - 1] == sequence2[i - 1] else DIFF)
                            # Value derived from DIAGONAL cell
                            ]
            scoringMatrix[i][j] = max(matrixValues)  # Take maximum value among three (Up, Left, Diagonal)

    return scoringMatrix

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

    """
    Read FASTA file and obtain both sequences
    Create output file if not exist
    """

    s1Fasta = open(args.sequence1, "r")
    s2Fasta = open(args.sequence2, "r")
    outputFile = open(args.outputFile, "w+")

    sequence1 = s1Fasta.read().replace("\n", "")
    sequence2 = s2Fasta.read().replace("\n", "")
    sequence1 = sequence1.replace(" ","")
    sequence2 = sequence2.replace(" ","")
    sequence1 = sequence1.upper()
    sequence2 = sequence2.upper()

    ##############################################################
    # Checking MAX_NUM_PATHS and Whether Sequences are Too Large #
    ##############################################################

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
        scoringMatrix = buildMatrix(sequence1, sequence2)
        print("Left Side Sequence: ", sequence2)
        print("Up side Sequence: ", sequence1)
        print("Maximum Paths: ", MAX_NUMBER_PATHS)
        print("\nScoring Matrix:\n", scoringMatrix)

        sequences = Alignment.needlemanWunschAlign(sequence1, sequence2, scoringMatrix)
        match, score, alignedSequence1, alignedSequence2 = finalize(sequences)

        """Writing final result in output file"""

        for i in range(0, len(match)):
            print("Alignment(Path) ", i+1, ":\n", alignedSequence1[i],"\n",match[i],"\n",alignedSequence2[i],"\nScore: ",score[i],"\n\n")
            outputFile.write("Alignment %s:\n%s\n%s\n%s\nScore: %s\n\n" % ((i + 1), alignedSequence1[i], match[i], alignedSequence2[i], score[i]))

        outputFile.close()

