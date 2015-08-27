# DISF test 1

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import DynamicIncoherentStructureFactor_serial

parameters = {}

parameters['version'] = "3.0.9"
parameters['estimate'] = "no"

parameters['weights'] = 'incoherent'
parameters['timeunits'] = 'ps'
parameters['frequencyunits'] = 'THz'
parameters['qunits'] = 'nm^-1'
parameters['output'] = os.path.join(tempfile.gettempdir(),"DISF_Current.nc")
parameters['pyroserver'] = 'monoprocessor'

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['timeinfo'] = '1:49:1'
parameters['subset'] = 'filename %s' % os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1_1.nms')
parameters['deuteration'] = None
parameters['resolution'] = 8.07276703
parameters['qvectors'] = {"qgeometry": "spatial", "qvectorspershell" : 500, "qshellwidth" : 1.0, "qshellvalues" : "3.:10.:1."}

analysis = DynamicIncoherentStructureFactor_serial(parameters)
analysis.runAnalysis()
