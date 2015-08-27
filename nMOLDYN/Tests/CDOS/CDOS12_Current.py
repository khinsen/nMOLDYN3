# CDOS test 12

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

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest2.nc')
parameters['timeinfo'] = "2:20:2"
parameters['weights'] = 'equal'
parameters['subset'] = 'object objectname P1960 misc backbone,sidechains'
parameters['deuteration'] = 'object objectname P1960 misc sidechains'
parameters['differentiation'] = 1
parameters['resolution'] = 10.76368937

analysis = CartesianDensityOfStates_serial(parameters)
analysis.runAnalysis()
