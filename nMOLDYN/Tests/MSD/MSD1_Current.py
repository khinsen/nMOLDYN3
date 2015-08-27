# MSD test 1

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
parameters['timeinfo'] = '1:49:1'
parameters['weights'] = 'mass'
parameters['subset'] = 'filename %s' % os.path.join(GVAR['nmoldyn_trajs'],'TrajectoryTest1_1.nms')
parameters['deuteration'] = None

analysis = MeanSquareDisplacement_serial(parameters)
analysis.runAnalysis()
