# A few (only a few I swear !!!) global variables.
  
import os
import platform

import MMTK

GVAR = {}

# The path to MMTK package.
GVAR['mmtk_path'] = MMTK.__path__[0]

# The path to nMOLDYN package.
GVAR['nmoldyn_path'] = os.path.split(__file__)[0]

# The path to nMOLDYN tests directory.
GVAR['nmoldyn_tests'] = os.path.join(GVAR['nmoldyn_path'], "Tests")

# The path to nMOLDYN analysis database.
GVAR['nmoldyn_analysis'] = os.path.join(GVAR['nmoldyn_path'], "Database", "Analysis")

# The path to nMOLDYN trajectories.
GVAR['nmoldyn_trajs'] = os.path.join(GVAR['nmoldyn_path'], "Trajectories")

# The path to pMOLDYN script.
GVAR['pmoldyn_path'] = os.path.join(GVAR['nmoldyn_tests'], "nMOLDYN_Reference", "pMoldyn.py")

# The path to the trajectory currently loaded.
GVAR['current_traj'] = None

# The list of the loaded trajectories.
GVAR['loaded_trajs'] = []

# The path to nMOLDYN jobs directory.
if platform.system().upper() == 'WINDOWS':
    GVAR['nmoldyn_jobs'] = os.path.join(os.environ['USERPROFILE'], 'Application Data', 'nMOLDYN', 'nmoldyn_jobs')
    
else:
    GVAR['nmoldyn_jobs'] = os.path.join(os.environ['HOME'], '.nmoldyn_jobs')

del os
del platform
del MMTK