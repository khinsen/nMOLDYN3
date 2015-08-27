# DISF test 9

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import DynamicIncoherentStructureFactor_serial

parameters = {}

parameters['version'] = "3.0.9"
parameters['estimate'] = "no"

# Constant parameters all over the tests.
parameters['weights'] = 'incoherent'
parameters['timeunits'] = 'ps'
parameters['frequencyunits'] = 'THz'
parameters['qunits'] = 'nm^-1'
parameters['output'] = os.path.join(tempfile.gettempdir(),"DISF_Current.nc")
parameters['pyroserver'] = 'monoprocessor'


parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')

parameters['timeinfo'] = '1:19:1'
parameters['subset'] = 'object objectname P892 misc sidechains'
parameters['deuteration'] = 'object objectname P892 atomelement hydrogen'
parameters['resolution'] = 86.10951495
parameters['qvectors'] = {"qgeometry": "spatial", "qvectorspershell" : 500, "qshellwidth" : 0.8, "qshellvalues" : "7.5:9.5:1.0"}

analysis = DynamicIncoherentStructureFactor_serial(parameters)
analysis.runAnalysis()

