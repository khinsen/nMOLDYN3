# ADOS test 10

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import AngularDensityOfStates_serial

parameters = {}

parameters['version'] = "3.0.9"
parameters['estimate'] = "no"

parameters['projection'] = None
parameters['output'] = os.path.join(tempfile.gettempdir(),"ADOS_Current.nc")
parameters['pyroserver'] = 'monoprocessor'
parameters['normalize'] = 'no'
parameters['frequencyunits'] = 'THz'
parameters['stepwiserbt'] = 'no'

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest3.nc')
parameters['timeinfo'] = "2:45:3"
parameters['group'] = 'object objectname Water atomelement * groupinglevel cluster'
parameters['referenceframe'] = 2
parameters['differentiation'] = 5
parameters['resolution'] = 18.45203892

analysis = AngularDensityOfStates_serial(parameters)
analysis.runAnalysis()
