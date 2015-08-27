import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

title = 'Elastic Incoherent Scattering Function'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
units_q = 1.0
weights = 'incoherent'
log_file = 'logfile.log'
output_files = {'eisf': os.path.join(tempfile.gettempdir(),'test.plot')}
q_vector_set = ([3.0, 5.0, 7.0, 9.0], 0.80000000000000004, 500, None)
deuter = None
atoms = {'Protein.0': ['Methyl']}
time_info = (0, 19, 2)
