import os
import sys
sys.path.insert(0, '/home/cs/pellegrini/nMOLDYN/development')

from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Core.IO import convertASCIIToNetCDF

pMoldyn='/home/cs/pellegrini/nMOLDYN/nMOLDYN2.2.5/nMoldyn/bin/pMoldyn'

if len(sys.argv) > 1:
    selectedTests = [int(v) for v in sys.argv[1:]]
else:
    selectedTests = range(1,13)

for i in selectedTests:
    
    LogMessage('info', "Run test %d" % i, ['console'])
    
    os.system(pMoldyn + ' --isf --input=DISF%s_Reference.inp' % i)
    os.system('mv DISF_Fqt_test.nc Fqt%s_Reference.nc' % i)
    os.system('mv DISF_Sqw_test.nc Sqw%s_Reference.nc' % i)
