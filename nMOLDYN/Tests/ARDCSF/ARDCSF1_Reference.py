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
q_vector_set = ([3.0, 4.0, 5.0], 1.0, 500, None)
deuter = None
ar_order = 5
atoms_pdb = os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1_1.pdb")
time_info = (0, 10, 1)
