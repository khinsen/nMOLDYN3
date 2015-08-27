# AVACF test 10

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import AngularVelocityAutoCorrelationFunction_serial

parameters = {}

parameters['version'] = "3.0.9"
parameters['estimate'] = "no"

parameters['timeunits'] = 'ps'
parameters['projection'] = None
parameters['output'] = os.path.join(tempfile.gettempdir(),"AVACF_Current.nc")
parameters['pyroserver'] = 'monoprocessor'
parameters['frequencyunits'] = 'THz'
parameters['stepwiserbt'] = 'no'

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest3.nc')
parameters['timeinfo'] = "2:45:3"
parameters['group'] = 'object objectname Water atomelement * groupinglevel cluster'
parameters['referenceframe'] = 2
parameters['differentiation'] = 5

analysis = AngularVelocityAutoCorrelationFunction_serial(parameters)
analysis.runAnalysis()
