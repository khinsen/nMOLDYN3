import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

units_length = 1.0
title = 'Mean-Square Displacement'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
log_file = 'logfile.log'
projection_vector = None
output_files = {'msd': os.path.join(tempfile.gettempdir(),'test.plot')}
deuter = None
atoms = {'Protein.0': ['Oxygen', 'Nitrogen']}
weights = 'none'
time_info = (1, 20, 2)
