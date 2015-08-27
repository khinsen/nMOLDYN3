# CDOS test 4

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
parameters['timeinfo'] = '17:39:4'
parameters['weights'] = 'equal'
parameters['subset'] = 'object objectname P892 chainname *'
parameters['deuteration'] = 'object objectname P892 chemfragment methyl'
parameters['resolution'] = 77.49856345

analysis = CartesianDensityOfStates_serial(parameters)
analysis.runAnalysis()
