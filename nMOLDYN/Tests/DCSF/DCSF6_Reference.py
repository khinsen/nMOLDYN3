import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

title = 'Coherent Scattering Function'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
units_q = 1.0
weights = 'coherent'
frequency_points = 1000000
log_file = 'logfile.log'
units_frequency = 1.0
output_files = {'csf': os.path.join(tempfile.gettempdir(),'DCSF_Reference.nc')}
q_vector_set = ([3.0, 4.0, 6.0, 7.0, 9.0], 1.0, 500, None)
deuter = {'Protein.0': ['Methyl']}
atoms = {'Protein.0': ['*']}
ft_window = 40.0
time_info = (4, 23, 2)
