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
q_vector_set = ([4.0, 5.0, 6.0, 7.0, 8.0, 9.0], 0.90000000000000002, 500, None)
deuter = None
atoms = {'Protein.0': ['Hydrogen']}
ft_window = 10.0
time_info = (5, 29, 3)
