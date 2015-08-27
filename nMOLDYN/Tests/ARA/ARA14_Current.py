# ARA test 14

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import AutoRegressiveAnalysis_serial

parameters = {}

parameters['output'] = os.path.join(tempfile.gettempdir(),"ARA_Current.nc")
parameters['pyroserver'] = 'monoprocessor'
parameters['projection'] = None

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest3.nc')
parameters['timeinfo'] = "2:45:3"
parameters['weights'] = 'equal'
parameters['subset'] = 'object objectname Water atomelement oxygen'
parameters['deuteration'] = None
parameters['differentiation'] = 1
parameters['armodelorder'] = 13

analysis = AutoRegressiveAnalysis_serial(parameters)
analysis.runAnalysis()
