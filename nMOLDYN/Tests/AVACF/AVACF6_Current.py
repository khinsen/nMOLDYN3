# AVACF test 6

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

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['timeinfo'] = "2:20:2"
parameters['group'] = 'object objectname P892 restype Val groupinglevel residue'
parameters['referenceframe'] = 2
parameters['differentiation'] = 4

analysis = AngularVelocityAutoCorrelationFunction_serial(parameters)
analysis.runAnalysis()
