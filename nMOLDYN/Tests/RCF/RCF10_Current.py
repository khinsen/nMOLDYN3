# RCF test 10

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
parameters['timeinfo'] = "2:33:2"
parameters['group'] = 'object objectname P892 restype Pro,Gln,Glu groupinglevel residue'
parameters['referenceframe'] = 2
parameters['wignerindexes'] = '0,0,0'

analysis = ReorientationalCorrelationFunction_serial(parameters)
analysis.runAnalysis()
