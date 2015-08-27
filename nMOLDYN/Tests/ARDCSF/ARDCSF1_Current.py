# ARDCSF test 1

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import AutoRegressiveDynamicCoherentStructureFactor_serial

parameters = {}

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['weights'] = 'coherent'
parameters['output'] = os.path.join(tempfile.gettempdir(),"ARDCSF_Current.nc")
parameters['pyroserver'] = 'monoprocessor'

parameters['timeinfo'] = '1:10:1'
parameters['subset'] = 'filename %s' % os.path.join(GVAR['nmoldyn_trajs'],'TrajectoryTest1_1.nms')
parameters['deuteration'] = None
parameters['armodelorder'] = 5
parameters['qvectors'] = {"qgeometry": "spatial", "qvectorspershell" : 500, "qshellwidth" : 1.0, "qshellvalues" : "3.:5.:1."}

analysis = AutoRegressiveDynamicCoherentStructureFactor_serial(parameters)
analysis.runAnalysis()
