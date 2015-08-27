# Run the complete test suite
#
# Written by Eric Pellegrini.
#

import unittest

import MSD_tests
import CVACF_tests
import CDOS_tests

def suite():
    test_suite = unittest.TestSuite()
    test_suite.addTests(MSD_tests.suite())
    test_suite.addTests(CVACF_tests.suite())
    test_suite.addTests(CDOS_tests.suite())
    return test_suite
    
if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
    
