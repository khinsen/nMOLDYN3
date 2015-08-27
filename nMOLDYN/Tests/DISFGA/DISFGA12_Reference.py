import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

title = 'Incoherent Scattering Function (Gaussian approx.)'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest3.nc")]
units_q = 1.0
weights = 'incoherent'
frequency_points = 1000000
log_file = 'logfile.log'
units_frequency = 1.0
output_files = {'isf': os.path.join(tempfile.gettempdir(),'DISFGA_Reference.nc')}
time_info = (0, 19, 1)
atoms = {'Water': ['Oxygen', 'Hydrogen']}
ft_window = 8.0
deuter = {'Water': ['Hydrogen']}
q_vector_set = ([5.5999999999999996, 6.4000000000000004, 7.2999999999999998, 10.5, 20.399999999999999], 0.59999999999999998, 500, None)
