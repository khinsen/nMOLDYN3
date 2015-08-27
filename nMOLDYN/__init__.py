#
# Package information
#


# Add shared library path to sys.path
import os
import sys

# Adds the path for the compiled modules to the system path.
sys.path.append(os.path.join(os.path.split(__file__)[0], sys.platform))

from nMOLDYN.Preferences import loadPreferencesFile, setPreferences

# The user configuration file is read.
prefs = loadPreferencesFile()
setPreferences(prefs)

del loadPreferencesFile
del os
del prefs
del setPreferences
del sys
