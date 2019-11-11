import unittest
from Alignment import *
from NeedlemanWunsch import *

class TestNeedleAlignment(unittest.TestCase):

    """Works when config file is set according to given value, otherwise make changes in example according to config file"""
    def setUp(self):
        self.sequence1 = "SUM"
        self.sequence2 = "SAM"
        self.GAP_PENALTY = -2
        self.SAME = -5
        self.DIFF = 5
        self.MAX_SEQ_LENGTH = 8
        self.MAX_NUMBER_PATHS = 1
        self.scoreMatrix = [[0, -2, -4, -6],
                            [-2, 5, 3, 1],
                            [-4, 3, 1, -1],
                            [-6, 1, -1, 6]]

    def testSequenceLength(self):
        self.assertLessEqual(len(self.sequence1), self.MAX_SEQ_LENGTH, "Sequence1 is too long.")
        self.assertLessEqual(len(self.sequence2), self.MAX_SEQ_LENGTH, "Sequence2 is too long")

    def testScoreMatrix(self):
        expectedMatrix = self.scoreMatrix

        resultMatrix = buildMatrix(self.sequence1, self.sequence2)

        for i in range(0, len(self.sequence2)):
            for j in range(0, len(self.sequence1)):
                self.assertEqual(resultMatrix[i][j],
                                 expectedMatrix[j][i],
                                 "Should return correct score matrix.")

    def testAlign(self):
        expectedAlign = [['S-UM', 'SU-M'], ['SA-M', 'S-AM']]
        actualAlign = needlemanWunschAlign(self.sequence1, self.sequence2, self.scoreMatrix)
        self.assertListEqual(actualAlign, expectedAlign, "Alignment is not proper")

    def testFinalScore(self):
        expectedScore = 6
        match, actualScore, finalAseq1, finalAseq2 = finalize([['S-UM', 'SU-M'], ['SA-M', 'S-AM']])
        self.assertEqual(actualScore[0], expectedScore, "Wrong Calculation and wrong Alignment")

    def testMaxPaths(self):
        match, actualScore, finalAseq1, finalAseq2 = finalize([['S-UM', 'SU-M'], ['SA-M', 'S-AM']])
        self.assertLessEqual(len(finalAseq1), self.MAX_NUMBER_PATHS, "Too many paths available")

    def testFinalPaths(self):
        match, actualScore, finalAseq1, finalAseq2 = finalize([['S-UM', 'SU-M'], ['SA-M', 'S-AM']])
        expectedPathSeq1 = ['S-UM']
        expectedPathSeq2 = ['SA-M']

        self.assertEqual(finalAseq1, expectedPathSeq1, "Wrong aligned sequence1")
        self.assertEqual(finalAseq2, expectedPathSeq2, "Wrong aligned sequence2")


if __name__ == '__main__':
    unittest.main()
