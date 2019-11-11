# Needleman-Wunsch
Global Alignment of FASTA sequences

USAGE:

python NeedlemanWunsch.py sequence1 sequence2 outputFile --gapPenalty negativeInteger --match positiveInteger --mismatch negativeInteger

where,
sequence1 --> path and name of the file for sequence1
sequence2 --> path and name of the file for sequence2
outputFile --> path and name for the output file

optional attributes:
--gapPenalty
--match
--mismatch

Example:
(for pyCharm:)
script path --> path of NeedlemanWunsch.py
parameters --> sequence1.fasta sequence2.fasta output.txt --gapPenalty -2 --match 5 --mismatch 5
python interpreter --> python 3.7


for tests:
1. script path --> path of TestAlignFunctions.py
   python interpreter --> python 3.7

2. script path --> path of TestNeedleWunschExample.py
   python interpreter --> python 3.7



# GLOBAL ALLIGNMENT ALGORITHM: 

    H(i,j) = max{ H(i-1, j-1) + S(A(i),B(j)), H(i-1, j) + GP, H(i, j-1) + GP }
    S(A(i),B(j)):
      if A(i) == B(j):-
        return match
      else:
        return mismatch

EXAMPLE: (seq1 = "SUM", seq2="SAM", GP = -2, match = 5, mismatch = 5)

          S  U  M
        0 -2 -4 -6
     S -2  5  3  1
     A -4  3  1 -1
     M -6  1 -1  6
     
   (SCORING MATRIX)
  
  
  Possible alignments:
  S - U M
  |     |
  S A - M
  
  
  S U - M
  |     |
  S - A M
  
  Score: 6
  
  
