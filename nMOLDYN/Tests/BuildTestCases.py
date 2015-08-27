# This python script creates automatically the test cases for a given analysis.
# It needs:
#   -the file 'testContents.txt' located in REFERENCES\ANALYSIS directory where
# ANALYSIS is the name of the analysis for which the test cases should be built (i.e. MSD, VACF ...).

import os
import sys

from nMOLDYN.GlobalVariables import GVAR

sys.path.insert(0, GVAR['nmoldyn_path'])

analysis = sys.argv[1]

testContents = open(os.path.join(GVAR['nmoldyn_path'], 'Tests', analysis, 'TestsContents.py'),'r')
exec testContents
testContents.close()

for t in range(len(test)):
    tInd = str(t + 1)
    referenceTestFile = open(os.path.join(GVAR['nmoldyn_path'],'Tests',analysis, analysis + tInd +'_Reference.inp'),'w')
    referenceTestFile.write('from MMTK import *\n')
    for key, val in template['REF'].items() + test[t]['REF'].items():
        if isinstance(val, str):
            referenceTestFile.write('%s = %s\n' % (key, repr(val)))
        else:
            referenceTestFile.write('%s = %s\n' % (key, val))
    referenceTestFile.close()
