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
q_vector_set = ([3.0, 4.0, 5.0, 6.0], 1.8, 500, None)
deuter = None
atoms = {'Protein.0': ['Sulfur']}
ft_window = 10.0
time_info = (3, 39, 4)
