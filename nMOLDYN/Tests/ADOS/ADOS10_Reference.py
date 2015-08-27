import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

title = 'Angular Velocity Autocorrelation Function'
ft_window = 20.0
time_steps = None
frequency_points = 1000000
log_file = 'logfile.log'
units_frequency = 1.0
projection_vector = None
output_files = {'savacf': os.path.join(tempfile.gettempdir(),'test.plot')}
differentiation = 'order 5'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest3.nc")]
reference = [{'Water': {'*': None}}]
groups = [{'Water': ['*']}]
time_info = (1, 45, 3)
