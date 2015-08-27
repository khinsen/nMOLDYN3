# DISFG test 11

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

parameters['trajectory'] = os.path.join(GVAR['nmoldyn_trajs'], 'TrajectoryTest3.nc')
parameters['timeinfo'] = '6:19:2'
parameters['subset'] = 'object objectname Water atomelement hydrogen'
parameters['deuteration'] = None
parameters['qshellvalues'] = "5.:12.:1.0"
parameters['resolution'] = 129.1610434

analysis = DynamicIncoherentStructureFactorGaussianApproximation_serial(parameters)
analysis.runAnalysis()
