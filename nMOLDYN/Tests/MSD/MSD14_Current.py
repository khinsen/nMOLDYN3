# MSD test 14

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import MeanSquareDisplacement_serial

parameters = {}

parameters['version'] = "3.0.9"
parameters['estimate'] = "no"

parameters['timeunits'] = 'ps'
parameters['distanceunits'] = 'nm'
parameters['projection'] = None
parameters['output'] = os.path.join(tempfile.gettempdir(),"MSD_Current.nc")
parameters['pyroserver'] = 'monoprocessor'

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['timeinfo'] = "2:40:2"
parameters['weights'] = 'equal'
parameters['subset'] = 'object objectname P892 chemfragment c_alphas'
parameters['deuteration'] = None

analysis = MeanSquareDisplacement_serial(parameters)
analysis.runAnalysis()
