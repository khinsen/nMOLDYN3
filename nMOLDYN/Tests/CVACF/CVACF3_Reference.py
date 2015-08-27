import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

units_length = 1.0
title = 'Velocity Auto-Correlation Function'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
log_file = 'logfile.log'
projection_vector = None
output_files = {'vacf': os.path.join(tempfile.gettempdir(),'test.plot')}
deuter = None
atoms = {'Protein.0': ['Methyl']}
weights = 'none'
time_info = (8, 29, 2)
