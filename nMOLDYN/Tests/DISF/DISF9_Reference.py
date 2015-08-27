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
q_vector_set = ([7.5, 8.5, 9.5], 0.80000000000000004, 500, None)
deuter = {'Protein.0': ['Hydrogen']}
atoms = {'Protein.0': ['SideChain']}
ft_window = 20.0
time_info = (0, 19, 1)
