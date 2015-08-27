import os
import sys
sys.path.insert(0, '/home/cs/pellegrini/nMOLDYN/development')

from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Core.IO import convertASCIIToNetCDF

pMoldyn='/home/cs/pellegrini/nMOLDYN/nMOLDYN2.2.5/nMoldyn/bin/pMoldyn'

if len(sys.argv) > 1:
    selectedTests = [int(v) for v in sys.argv[1:]]
else:
    selectedTests = range(1,15)

a = 'ARA'

for i in selectedTests:
    LogMessage('info', "Run test %d" % i, ['console'])
    if i <= 7:
        os.system(pMoldyn + ' --ar-vel --input=%s%s_Reference.inp' % (a,i))
    else:
        os.system(pMoldyn + ' --ar-xyz --input=%s%s_Reference.inp' % (a,i))

    os.system('mv testDOS.plot DOS%s_Reference.plot' % i)
    os.system('mv testMSD.plot MSD%s_Reference.plot' % i)
    os.system('mv testVACF.plot VACF%s_Reference.plot' % i)
    os.system('mv testPARAMETERS.plot PARAMETERS%s_Reference.plot' % i)
    os.system('mv testMEMORY.plot MEMORY%s_Reference.plot' % i)
    convertASCIIToNetCDF('DOS%s_Reference.plot' % i, 'DOS%s_Reference.nc' % i)
    convertASCIIToNetCDF('MSD%s_Reference.plot' % i, 'MSD%s_Reference.nc' % i)
    convertASCIIToNetCDF('VACF%s_Reference.plot' % i, 'VACF%s_Reference.nc' % i)
    convertASCIIToNetCDF('MEMORY%s_Reference.plot' % i, 'MEMORY%s_Reference.nc' % i)
