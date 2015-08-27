import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

title = 'Angular Velocity Autocorrelation Function'
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
time_steps = None
log_file = 'logfile.log'
projection_vector = None
output_files = {'avacf': os.path.join(tempfile.gettempdir(),'test.plot')}
differentiation = 'order 4'
reference = [{'Protein.0 SideChain': {'leucine': None, 'methionine': None}}]
groups = [{'Protein.0 SideChain': ['leucine', 'methionine']}]
time_info = (1, 40, 2)
