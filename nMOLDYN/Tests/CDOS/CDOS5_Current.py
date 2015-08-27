# CDOS test 5

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import CartesianDensityOfStates_serial

parameters = {}

parameters['version'] = "3.0.9"
parameters['estimate'] = "no"

parameters['timeunits'] = 'ps'
parameters['projection'] = None
parameters['output'] = os.path.join(tempfile.gettempdir(),"CDOS_Current.nc")
parameters['pyroserver'] = 'monoprocessor'
parameters['normalize'] = 'no'
parameters['differentiation'] = 1
parameters['frequencyunits'] = 'THz'

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['timeinfo'] = '2:20:2'
parameters['weights'] = 'equal'
parameters['subset'] = 'object objectname P892 atomelement oxygen,nitrogen'
parameters['deuteration'] = None
parameters['resolution'] = 430.54757474

analysis = CartesianDensityOfStates_serial(parameters)
analysis.runAnalysis()
