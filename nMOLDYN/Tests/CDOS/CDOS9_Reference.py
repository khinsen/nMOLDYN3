import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

units_length = 1.0
title = 'Density Of States'
frequency_points = 1000000
log_file = 'logfile.log'
units_frequency = 1.0
projection_vector = None
output_files = {'dos': os.path.join(tempfile.gettempdir(),'test.plot')}
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest2.nc")]
time_info = (1, 20, 1)
atoms = {'Protein.0': ['Hydrogen']}
ft_window = 4.0
deuter = None
weights = 'none'
