# ARA test 13

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import AutoRegressiveAnalysis_serial

parameters = {}

parameters['output'] = os.path.join(tempfile.gettempdir(),"ARA_Current.nc")
parameters['pyroserver'] = 'monoprocessor'
parameters['projection'] = None

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest3.nc')
parameters['timeinfo'] = "4:20:2"
parameters['weights'] = 'equal'
parameters['subset'] = 'object objectname Water atomelement oxygen,hydrogen'
parameters['deuteration'] = 'object objectname Water atomelement hydrogen'
parameters['differentiation'] = 1
parameters['armodelorder'] = 4

analysis = AutoRegressiveAnalysis_serial(parameters)
analysis.runAnalysis()
