import unittest
from NeedlemanWunsch import *

class TestNeedleAlignment(unittest.TestCase):
    def testExample1(self):
        if MAX_NUMBER_PATHS == 2 and GAP_PENALTY == -2 and SAME == 5 and DIFF == -5:
            """ Given """
            sequence1 = "GCT"
            sequence2 = "GAAT"

            """ When"""
            actual_matrix, actual_alignment_seq1, actual_alignment_seq2, match = align_sequences(sequence1, sequence2)
            actual_score = actual_matrix[len(sequence1)][len(sequence2)]
            expected_matrix = [[0, -2, -4, -6, -8],
                               [-2, 5, 3, 1, -1],
                               [-4, 3, 1, -1, -3],
                               [-6, 1, -1, -3, 4]]
            expected_score = 4
            expected_alignment_seq1 = ["GC--T", "G-C-T"]
            expected_alignment_seq2 = ["G-AAT", "GA-AT"]

            """ Then """
            for i in range(0, len(sequence1)):
                for j in range(0, len(sequence2)):
                    self.assertEqual(actual_matrix[i][j],
                                     expected_matrix[i][j],
                                     "Should return correct score matrix.")

            self.assertEqual(actual_score, expected_score, "Wrong Score")
            self.assertListEqual(actual_alignment_seq1, expected_alignment_seq1, "Alignment1 is not proper")
            self.assertListEqual(actual_alignment_seq2, expected_alignment_seq2, "Alignment2 is not proper")
        else:
            stderr.write("\nSet config variables MAX_NUMBER_PATHS = 2, GAP_PENALTY = -2, SAME = 5, DIFF = -5\n\n")
            exit(1)


if __name__ == '__main__':
    unittest.main()
