import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

units_length = 1.0
title = 'Autoregressive Model'
ar_precision = None
time_steps = None
frequency_points = None
log_file = 'logfile.log'
units_frequency = 1.0
projection_vector = None
output_files = {'ara': os.path.join(tempfile.gettempdir(),'test.plot')}
trajectory = [os.path.join(GVAR['nmoldyn_trajs'],'TrajectoryTest3.nc')]
time_info = (3, 20, 2)
ar_order = 4
atoms = {'Water': ['Oxygen', 'Hydrogen']}
weights = 'none'
deuter = {'Water': ['Hydrogen']}
