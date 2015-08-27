# ARA test 1

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import AutoRegressiveAnalysis_serial

parameters = {}

parameters['output'] = os.path.join(tempfile.gettempdir(),"ARA_Current.nc")
parameters['pyroserver'] = 'monoprocessor'
parameters['projection'] = None

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['timeinfo'] = '1:19:1'
parameters['weights'] = 'mass'
parameters['subset'] = 'filename %s' % os.path.join(GVAR['nmoldyn_trajs'],'TrajectoryTest1_2.nms')
parameters['deuteration'] = None
parameters['differentiation'] = 1
parameters['armodelorder'] = 10

analysis = AutoRegressiveAnalysis_serial(parameters)
analysis.runAnalysis()
