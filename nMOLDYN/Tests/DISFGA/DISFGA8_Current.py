# DISFG test 8

import os
import tempfile

from nMOLDYN.GlobalVariables import GVAR

from nMOLDYN.Analysis.Templates import DynamicIncoherentStructureFactorGaussianApproximation_serial

parameters = {}

parameters['version'] = "3.0.9"
parameters['estimate'] = "no"

parameters['weights'] = 'incoherent'
parameters['timeunits'] = 'ps'
parameters['frequencyunits'] = 'THz'
parameters['qunits'] = 'nm^-1'
parameters['output'] = os.path.join(tempfile.gettempdir(),"DISFGA_Current.nc")
parameters['pyroserver'] = 'monoprocessor'

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest1.nc')
parameters['timeinfo'] = '4:39:4'
parameters['subset'] = 'object objectname P892 atomelement sulfur'
parameters['deuteration'] = None
parameters['qshellvalues'] = "3.:6.:1.0"
parameters['resolution'] = 96.87320432

analysis = DynamicIncoherentStructureFactorGaussianApproximation_serial(parameters)
analysis.runAnalysis()

