# DCSF test 12

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
parameters['timeinfo'] = '1:19:1'
parameters['subset'] = 'object objectname P892 atomelement carbon,hydrogen,nitrogen,oxygen,sulfur'
parameters['deuteration'] = None
parameters['resolution'] = 215.27378737
parameters['qvectors'] = {"qgeometry": "spatial", "qvectorspershell" : 500, "qshellwidth" : 0.6, "qshellvalues" : (3.2, 4.2, 5.6, 6.4, 7.3)}

analysis = DynamicCoherentStructureFactor_serial(parameters)
analysis.runAnalysis()


