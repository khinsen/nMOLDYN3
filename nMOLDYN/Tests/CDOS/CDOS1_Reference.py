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
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
time_info = (0, 49, 1)
ft_window = 20.0
deuter = None
weights = 'mass'
atoms_pdb = os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1_2.pdb")
