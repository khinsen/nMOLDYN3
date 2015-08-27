import os
import sys
sys.path.insert(0, '/home/cs/pellegrini/nMOLDYN/development')

from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Core.IO import convertASCIIToNetCDF

pMoldyn='/home/cs/pellegrini/nMOLDYN/nMOLDYN2.2.5/nMoldyn/bin/pMoldyn'

if len(sys.argv) > 1:
    selectedTests = [int(v) for v in sys.argv[1:]]
else:
    selectedTests = range(1,11)

a = 'ADOS'

for i in selectedTests:
    LogMessage('info', "Run test %d" % i, ['console'])
    os.system(pMoldyn + ' --savacf --input=%s%s_Reference.inp' % (a,i))

    os.system('mv test.plot %s%s_Reference.plot' % (a,i))
    convertASCIIToNetCDF('%s%s_Reference.plot' % (a,i), '%s%s_Reference.nc' % (a,i))
