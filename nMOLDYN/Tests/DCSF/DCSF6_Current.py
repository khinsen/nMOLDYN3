# DCSF test 6

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import DynamicCoherentStructureFactor_serial

parameters = {}

parameters['version'] = "3.0.9"
parameters['estimate'] = "no"

parameters['weights'] = 'coherent'
parameters['timeunits'] = 'ps'
parameters['frequencyunits'] = 'THz'
parameters['qunits'] = 'nm^-1'
parameters['output'] = os.path.join(tempfile.gettempdir(),"DCSF_Current.nc")
parameters['pyroserver'] = 'monoprocessor'

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['timeinfo'] = '5:23:2'
parameters['subset'] = None
parameters['deuteration'] = 'object objectname P892 chemfragment methyl'
parameters['resolution'] = 43.05475747
parameters['qvectors'] = {"qgeometry": "spatial", "qvectorspershell" : 500, "qshellwidth" : 1.0, "qshellvalues" : [3.0, 4.0, 6.0, 7.0, 9.0]}

analysis = DynamicCoherentStructureFactor_serial(parameters)
analysis.runAnalysis()
