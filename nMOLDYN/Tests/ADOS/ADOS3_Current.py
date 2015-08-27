# ADOS test 3

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

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['timeinfo'] = '9:29:2'
parameters['group'] = 'object objectname P892 chemfragment methyl groupinglevel methyl'
parameters['referenceframe'] = 9
parameters['differentiation'] = 2
parameters['resolution'] = 61.99885076

analysis = AngularDensityOfStates_serial(parameters)
analysis.runAnalysis()
