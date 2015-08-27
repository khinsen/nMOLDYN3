# DISFG test 7

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
parameters['timeinfo'] = '1:19:1'
parameters['subset'] = 'object objectname P892 atomelement carbon'
parameters['deuteration'] = None
parameters['qshellvalues'] = "4.:8.:1.0"
parameters['resolution'] = 21.52737874

analysis = DynamicIncoherentStructureFactorGaussianApproximation_serial(parameters)
analysis.runAnalysis()
