import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

title = 'Incoherent Scattering Function'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
units_q = 1.0
weights = 'incoherent'
frequency_points = 1000000
log_file = 'logfile.log'
units_frequency = 1.0
output_files = {'isf': os.path.join(tempfile.gettempdir(),'DISF_Reference.nc')}
q_vector_set = ([3.2000000000000002, 4.2000000000000002, 5.5999999999999996, 6.4000000000000004, 7.2999999999999998], 0.59999999999999998, 500, None)
deuter = None
atoms = {'Protein.0': ['Oxygen', 'Nitrogen', 'Carbon', 'Sulfur', 'Hydrogen']}
ft_window = 8.0
time_info = (0, 19, 1)
