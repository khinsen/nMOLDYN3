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
    
    os.system(pMoldyn + ' --eisf --input=EISF%s_Reference.inp' % i)

    os.system('mv test.plot EISF%s_Reference.plot' % i)
    convertASCIIToNetCDF('EISF%s_Reference.plot' % i, 'EISF%s_Reference.nc' % i)
