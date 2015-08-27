import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

output_files = {'rcf': os.path.join(tempfile.gettempdir(),'test.plot')}
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
time_steps = None
log_file = 'logfile.log'
title = 'Reorientational Correlation Function'
rotation_coefficients = (8, 5, 4)
reference = [{'Protein.0 SideChain': {'alanine': None}}]
groups = [{'Protein.0 SideChain': ['alanine']}]
time_info = (1, 28, 2)
