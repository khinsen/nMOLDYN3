import subprocess
import sys
import os

removeProg = os.path.join(sys.prefix,'RemovenMOLDYN.exe')
logFile = os.path.join(sys.prefix, 'nMOLDYN-wininst.log')

subprocess.call([removeProg, '-u', logFile])
