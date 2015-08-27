# ARA test 6

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import AutoRegressiveAnalysis_serial

parameters = {}

parameters['output'] = os.path.join(tempfile.gettempdir(),"ARA_Current.nc")
parameters['pyroserver'] = 'monoprocessor'
parameters['projection'] = None

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['timeinfo'] = "2:20:2"
parameters['weights'] = 'mass'
parameters['subset'] = 'object objectname P892 atomelement nitrogen'
parameters['deuteration'] = None
parameters['differentiation'] = 1
parameters['armodelorder'] = 7

analysis = AutoRegressiveAnalysis_serial(parameters)
analysis.runAnalysis()
