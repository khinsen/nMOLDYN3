# RCF test 14

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import ReorientationalCorrelationFunction_serial

parameters = {}

parameters['version'] = "3.0.9"
parameters['estimate'] = "no"

parameters['timeunits'] = 'ps'
parameters['output'] = os.path.join(tempfile.gettempdir(),"RCF_Current.nc")
parameters['pyroserver'] = 'monoprocessor'
parameters['stepwiserbt'] = 'no'

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest3.nc')
parameters['timeinfo'] = "2:45:3"
parameters['group'] = 'object objectname Water atomelement * groupinglevel cluster'
parameters['referenceframe'] = 2
parameters['wignerindexes'] = '3,2,1'

analysis = ReorientationalCorrelationFunction_serial(parameters)
analysis.runAnalysis()
