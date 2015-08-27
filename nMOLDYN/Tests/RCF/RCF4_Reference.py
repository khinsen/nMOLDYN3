import os
import tempfile
from nMOLDYN.GlobalVariables import GVAR

output_files = {'rcf': os.path.join(tempfile.gettempdir(),'test.plot')}
trajectory = [os.path.join(GVAR['nmoldyn_trajs'], "TrajectoryTest1.nc")]
time_steps = None
log_file = 'logfile.log'
title = 'Reorientational Correlation Function'
rotation_coefficients = (5, 0, 0)
reference = [{'Protein.0 SideChain': {'alanine': None, 'threonine': None, 'glycine': None, 'valine': None, 'arginine': None, 'asparagine': None, 'glutamic_acid': None, 'leucine': None, 'methionine': None, 'glutamine': None, 'isoleucine': None, 'phenylalanine': None, 'proline': None, 'serine': None, 'tyrosine': None}}]
groups = [{'Protein.0 SideChain': ['alanine', 'glycine', 'proline', 'glutamine', 'arginine', 'glutamic_acid', 'leucine', 'serine', 'valine', 'methionine', 'isoleucine', 'asparagine', 'phenylalanine', 'threonine', 'tyrosine']}]
time_info = (16, 39, 4)
