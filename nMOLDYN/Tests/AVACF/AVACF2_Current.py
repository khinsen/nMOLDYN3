# AVACF test 2

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
parameters['timeinfo'] = '1:49:4'
parameters['group'] = 'object objectname P892 chemfragment methyl groupinglevel methyl'
parameters['referenceframe'] = 1
parameters['differentiation'] = 1

analysis = AngularVelocityAutoCorrelationFunction_serial(parameters)
analysis.runAnalysis()


