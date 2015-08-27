import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

output_files = {'rcf': os.path.join(tempfile.gettempdir(),'test.plot')}
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
time_steps = None
log_file = 'logfile.log'
title = 'Reorientational Correlation Function'
rotation_coefficients = (2, 1, 1)
reference = [{'Protein.0 SideChain': {'leucine': None, 'methionine': None}}]
groups = [{'Protein.0 SideChain': ['leucine', 'methionine']}]
time_info = (1, 40, 2)
