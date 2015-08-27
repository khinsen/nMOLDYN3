import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

output_files = {'rcf': os.path.join(tempfile.gettempdir(),'test.plot')}
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
time_steps = None
log_file = 'logfile.log'
title = 'Reorientational Correlation Function'
rotation_coefficients = (0, 0, 0)
reference = [{'Protein.0 BackBone': {'*': None}}]
groups = [{'Protein.0 BackBone': ['*']}]
time_info = (0, 19, 1)
