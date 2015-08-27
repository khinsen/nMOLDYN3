# DCSF test 5

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
parameters['timeinfo'] = '6:29:3'
parameters['subset'] = 'object objectname P892 atomelement hydrogen'
parameters['deuteration'] = None
parameters['resolution'] = 147.61631134
parameters['qvectors'] = {"qgeometry": "spatial", "qvectorspershell" : 500, "qshellwidth" : 0.9, "qshellvalues" : "4.:9.:1.0"}

analysis = DynamicCoherentStructureFactor_serial(parameters)
analysis.runAnalysis()
