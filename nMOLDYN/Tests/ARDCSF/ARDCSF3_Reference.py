import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

title = 'Coherent Scattering AR analysis'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
ar_precision = None
units_q = 1.0
weights = 'coherent'
frequency_points = 1000000
log_file = 'logfile.log'
units_frequency = 1.0
output_files = {'csf': os.path.join(tempfile.gettempdir(),'ARDCSF_Reference.nc')}
q_vector_set = ([3.0, 5.0, 7.0, 9.0], 0.80000000000000004, 500, None)
deuter = None
atoms = {'Protein.0': ['Methyl']}
ar_order = 3
time_info = (0, 13, 2)
