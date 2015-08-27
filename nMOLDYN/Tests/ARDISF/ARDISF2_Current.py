# ARDISF test 2

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import AutoRegressiveDynamicIncoherentStructureFactor_serial

parameters = {}

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['weights'] = 'incoherent'
parameters['output'] = os.path.join(tempfile.gettempdir(),"ARDISF_Current.nc")
parameters['pyroserver'] = 'monoprocessor'

parameters['timeinfo'] = '1:8:1'
parameters['subset'] = 'object objectname P892 atomelement sulfur'
parameters['deuteration'] = None
parameters['armodelorder'] = 4
parameters['qvectors'] = {"qgeometry": "spatial", "qvectorspershell" : 500, "qshellwidth" : 0.8, "qshellvalues" : "6.:8.:1.0"}

analysis = AutoRegressiveDynamicIncoherentStructureFactor_serial(parameters)
analysis.runAnalysis()
