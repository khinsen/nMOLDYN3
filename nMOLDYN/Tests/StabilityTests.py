"""
Contains the stability benchmarks for the common analysis between nMOLDYN 2 (2.2.5) and 3.
Each benchmark is made of several small tests corresponding to various input values. 
The output of these tests are:
    * stored in NetCDF output file for nMOLDYN 2
    * generated on the fly for nMOLDYN 3
A difference of a maximum of 1.0E-6 is tolerated for the difference between the results obtained with nMOLDYN 2 and 3.
"""

# The python distribution modules
import os
import subprocess
import sys
import tempfile
from timeit import default_timer
import unittest
from unittest import TestResult

# The ScientificPython modules
from Scientific import N
from Scientific.IO.NetCDF import NetCDFFile

# The nMOLDYN modules
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Core.Logger import LogMessage

TEMP_DIR = tempfile.gettempdir()
                                       
class StabilityTests(unittest.TestCase):
    
    def __init__(self, shortName, longName, pMoldynArg, id):
        
        unittest.TestCase.__init__(self)
        
        self.shortName = shortName

        self.longName = longName
        
        self.pMoldynArg = pMoldynArg
        
        self.id = id
        
        self.errors = {}
        
    def setUp(self):
        """Overides the TestCase.setUp method. Initialization of the test variables.
        """
                
        # The path specific to the selected benchmark.
        self.testPath = os.path.join(GVAR['nmoldyn_tests'], self.shortName)
        
    def runTest(self):
                
        # Some information about the test to run are displayed in the console and file loggers.
        LogMessage('info', 'ANALYSIS: %s --- TEST No: %s' % (self.longName, self.id),['console', 'file'])

        subprocess.call([sys.executable, GVAR['pmoldyn_path'], self.pMoldynArg, "--input", os.path.join(self.testPath, self.shortName + "%s_Reference.py" % self.id)])
        subprocess.call([sys.executable, os.path.join(self.testPath, self.shortName + "%d_Current.py" % self.id)])
        
        refFileName = os.path.join(TEMP_DIR, self.shortName + "_Reference.nc")
        curFileName = os.path.join(TEMP_DIR, self.shortName + "_Current.nc")

        refFile = NetCDFFile(refFileName, 'r')
        curFile = NetCDFFile(curFileName, 'r')
                                                
        for key, val in refFile.variables.items():
                                        
            if val.typecode() != 'c':
                                
                if not curFile.variables.has_key(key):
                    continue
                    
                refValue = val.getValue()
                curValue = curFile.variables[key].getValue()
                                                                                                                                                                                                                                                                                                    
                # Their difference is computed.
                errorMax = max(abs(N.ravel(refValue - curValue)))
                    
                self.errors[key] = errorMax
                                                                                                                                        
                # Throws an assertion error if the difference is bigger than 1.0E-6.
                self.assertAlmostEqual(errorMax, 0.0, 6)

        curFile.close()
        
        refFile.close()
        
        try:
            os.remove(refFileName)
            os.remove(curFileName)
        except:
            pass
            
class StabilityTestsResults(TestResult):

    def __init__(self):
        
        TestResult.__init__(self)
        
        self.duration = 0.0
        
        self.infos = {}
        
        self.failures = []
        self.errors = []
        self.successes = []

    def addFailure(self, test, err):
        
        self.failures.append(test.id)

        self.infos[test.id] = {"status" : "failure", "errors" : test.errors}

    def addError(self, test, err):
        
        self.errors.append(test.id)

        self.infos[test.id] = {"status" : "error", "errors" : test.errors}

    def addSuccess(self, test):
        
        self.successes.append(test.id)

        self.infos[test.id] = {"status" : "success", "errors" : test.errors}
        
    def summary(self):
                
        summ = ""
        summ += "%d tests run in %s seconds\n\n" % (self.testsRun, self.duration)
                
        summ += "Successfull tests: %s\n" % " ".join([str(v) for v in self.successes])
        summ += "Failing tests: %s\n" % " ".join([str(v) for v in self.failures])
        summ += "Errors: %s\n\n" % " ".join([str(v) for v in self.errors])
        
        summ += "Details:\n\n"
        
        for k in sorted(self.infos.keys()):

            summ += "Test %d:\n" % k
                        
            for kk in sorted(self.infos[k]["errors"].keys()):
                
                summ += "\t%s: %s\n" % (kk, self.infos[k]["errors"][kk])
                
            summ += "\n"
            
        return summ
                                                
if __name__ == '__main__':
    
    selectedTests = [n.upper() for n in sys.argv[1:]]
    
    pass
            

