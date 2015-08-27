# RCF test 12

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

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['timeinfo'] = "2:40:2"
parameters['group'] = 'object objectname P892 restype Leu,Met groupinglevel residue'
parameters['referenceframe'] = 2
parameters['wignerindexes'] = '2,1,1'

analysis = ReorientationalCorrelationFunction_serial(parameters)
analysis.runAnalysis()
