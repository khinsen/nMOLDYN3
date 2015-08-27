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
time_info = (5, 19, 2)
atoms = {'Water': ['Hydrogen']}
ft_window = 10.0
deuter = None
q_vector_set = ([5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0], 1.0, 500, None)
