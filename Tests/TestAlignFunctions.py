import unittest
from Alignment import *

class TestNeedleAlignment(unittest.TestCase):
    def setUp(self):
        self.s1 = 'H'
        self.s2 = 'H'
        self.s3 = 'I'
        self.row = 1
        self.col = 1

    def testDiaNoGap(self):
        expectedS1, expectedS2 = 'H', 'H'
        row, col, aS1, aS2 = dia(self.row, self.col,"", "", self.s1, self.s2)
        self.assertEqual(aS1, expectedS1, "Mistake in Sequence1 alignment while aligning element in diagonal way")
        self.assertEqual(aS2, expectedS2, "Mistake in Sequence2 alignment while aligning element in diagonal way")

    def testLeftGap(self):
        expectedS1 = 'H'
        expectedS2 = '-'
        row, col, aS1, aS2 = left(self.row, self.col,"", "", self.s1, self.s3)
        self.assertEqual(aS1, expectedS1, "Mistake in Sequence1 alignment while aligning element in left direction")
        self.assertEqual(aS2, expectedS2, "Mistake in Sequence1 alignment while aligning element in left direction")

    def testUpGap(self):
        expectedS1 = '-'
        expectedS2 = 'I'
        row, col, aS1, aS2 = up(self.row, self.col,"", "", self.s1, self.s3)
        self.assertEqual(aS1, expectedS1, "Mistake in Sequence1 alignment while aligning element in up direction")
        self.assertEqual(aS2, expectedS2, "Mistake in Sequence1 alignment while aligning element in up direction")

    def testAlignLastSymbol(self):
        expectedS1 = 'H-'
        expectedS2 = '-I'
        aS1, aS2, row, col = align("", "", self.row, self.col, self.s1, self.s3)
        self.assertEqual(aS1, expectedS1, "Mistake in putting gap in sequence1 after finish tracing up to the top most cell")
        self.assertEqual(aS2, expectedS2, "Mistake in putting gap in sequence2 after finish tracing up to the top most cell")

if __name__ == '__main__':
    unittest.main()

