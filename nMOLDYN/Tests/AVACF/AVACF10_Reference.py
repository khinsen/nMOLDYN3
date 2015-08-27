import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

title = 'Angular Velocity Autocorrelation Function'
time_steps = None
log_file = 'logfile.log'
projection_vector = None
output_files = {'avacf': os.path.join(tempfile.gettempdir(),'test.plot')}
differentiation = 'order 5'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest3.nc")]
reference = [{'Water': {'*': None}}]
groups = [{'Water': ['*']}]
time_info = (1, 45, 3)
