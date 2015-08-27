import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

title = 'Incoherent Scattering AR analysis'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
ar_precision = None
units_q = 1.0
weights = 'incoherent'
frequency_points = 1000000
log_file = 'logfile.log'
units_frequency = 1.0
output_files = {'isf': os.path.join(tempfile.gettempdir(),'ARDISF_Reference.nc')}
q_vector_set = ([3.0, 4.0, 5.0], 1.0, 500, None)
deuter = None
ar_order = 5
atoms_pdb = os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1_1.pdb")
time_info = (0, 10, 1)
