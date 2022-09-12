"""
Test cases for analysis modules.
    - testCorrelation1: test for the autocorrelation of a single time serie.
    - testCorrelation2: test for the autocorrelation of a two dimensional time serie.
    - testCorrelation3: another test for the autocorrelation of a two dimensional time serie.
"""

import unittest

from Scientific.N import array, zeros

from nMOLDYN.Mathematics.Analysis import correlation

class AnalysisTest(unittest.TestCase):

    """
    Tests Analysis module functionnalities.
    """

    def setUp(self):
        # The time serie for test 1
        self.series1 = array([0,1,0,4,6,5])
        # The theoretical result for test 1
        self.results1 = array([13.0,10.8,6.0,2.0,2.5,0.0])

        # The time serie for test 2
        self.series2 = array([
            [0,0,0,0,0,0,0,0,0,0],
            [1,1,1,1,1,1,1,1,1,1],
            [0,0,0,0,0,0,0,0,0,0],
            [4,4,4,4,4,4,4,4,4,4],
            [6,6,6,6,6,6,6,6,6,6],
            [5,5,5,5,5,5,5,5,5,5]
        ])
        # The theoretical result for test 2
        self.results2 = array([130.0,108.0,60.0,20.0,25.0,0.0])

        # The time serie for test 3
        self.series3 = array([
            [0,2,0],
            [1,1,-1],
            [0,1,0],
            [4,0,-4],
            [6,5,0],
            [5,0,3]
        ])

        # The theoretical result for test 3
        self.results3 = array([22.5,11.4,5.75,11.0/3,6,0])

    def testCorrelation1(self):
        corr = correlation(self.series1)
        errorMax = max(abs(corr-self.results1))
        self.assertAlmostEqual(errorMax, 0.0, 10)

    def testCorrelation2(self):
        corr = correlation(self.series2)
        errorMax = max(abs(corr-self.results2))
        self.assertAlmostEqual(errorMax, 0.0, 10)

    def testCorrelation3(self):
        corr = correlation(self.series3)
        errorMax = max(abs(corr-self.results3))
        self.assertAlmostEqual(errorMax, 0.0, 10)

if __name__ == '__main__':
    unittest.main()
