import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

output_files = {'rcf': os.path.join(tempfile.gettempdir(),'test.plot')}
time_steps = None
log_file = 'logfile.log'
title = 'Reorientational Correlation Function'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest3.nc")]
rotation_coefficients = (3, 2, 1)
reference = [{'Water': {'*': None}}]
groups = [{'Water': ['*']}]
time_info = (1, 45, 3)
