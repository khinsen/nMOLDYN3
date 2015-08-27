#!python

from glob import glob
import os
import platform
import re
import sys

from distutils.core import setup, Extension
from distutils.version import LooseVersion
import distutils.sysconfig
sysconfig = distutils.sysconfig.get_config_vars()
python_site_packages = distutils.sysconfig.get_python_lib()

def cmpver(a, b, pad = None):
    
    def fixup(i):
        try:
            return int(i)
        except ValueError:
            return i

    a = [fixup(i) for i in re.findall("\d+|\w+", str(a))]
    b = [fixup(i) for i in re.findall("\d+|\w+", str(b))]

    if pad is not None:
        if len(a) < len(b):
            [a.append(fixup(pad)) for i in range(len(a)+1,len(b)+1)]
        elif len(a) > len(b):
            [b.append(fixup(pad)) for i in range(len(b)+1,len(a)+1)]
            
    # -1 if a<b, 0 if a=b, 1 if a>b
    return cmp(a, b)

class Dummy(object):
    pass

pkginfo = Dummy()
execfile('nMOLDYN/__pkginfo__.py', pkginfo.__dict__)                               

# Check that we have python between 2.5.0 and 2.7.0 (inclusive).
try:
    python_minversion = LooseVersion(vstring = "2.5.0")
    python_maxversion = LooseVersion(vstring = "3.0.0")
    python_version = LooseVersion(vstring = '.'.join([str(v) for v in sys.version_info[0:3]]))
    
    if cmpver(python_version, python_minversion, pad = 0) == -1: raise 
    if cmpver(python_version, python_maxversion, pad = 0) ==  1: raise 
    
except :
    print 'nMOLDYN needs Python between %s and %s. Installation stopped.' % (python_minversion,python_maxversion)
    raise SystemExit

# Check that we have the minimal numpy version.
try:
    import numpy
    numpy_minversion = LooseVersion(vstring = "1.1.0")
    numpy_version = LooseVersion(vstring = numpy.__version__)
    
    if cmpver(numpy_version, numpy_minversion, pad = 0) == -1: raise
    
except :
    print 'nMOLDYN needs NumPy %s or higher. Installation stopped.' % numpy_minversion
    raise SystemExit

# Check that we have the minimal matplotlib version.
try:
    from matplotlib import __version__ as matplotlib_version
    matplotlib_minversion = LooseVersion(vstring = '0.98.0')
    matplotlib_version = LooseVersion(vstring = matplotlib_version)
    
    if cmpver(matplotlib_version, matplotlib_minversion, pad = 0) == -1: raise 

except :
    print 'nMOLDYN needs matplotlib %s or higher. The plotting functionnalities will be disabled.' % matplotlib_minversion

# Check that we have the minimal PyRO version.
try:
    from Pyro.constants import VERSION as pyro_version
    pyro_minversion = LooseVersion(vstring = '3.9.0')
    pyro_version = LooseVersion(vstring = pyro_version)
    
    if cmpver(pyro_version, pyro_minversion, pad = 0) == -1: raise 

except :
    print 'nMOLDYN needs pyro %s or higher. The parallel computing functionnalities will be disabled.' % pyro_minversion

# Check that we have the minimal Scientific version.
try:
    from Scientific import __version__ as scientific_version
    scientific_minversion = LooseVersion(vstring = '2.8.0')
    scientific_version = LooseVersion(vstring = scientific_version)
    
    if cmpver(scientific_version, scientific_minversion, pad = 0) == -1: raise 

except :
    print 'nMOLDYN needs ScientificPython %s or higher. Installation stopped.' % scientific_minversion
    raise SystemExit

# Check that we have the minimal MMTK version.
try:
    from MMTK import __version__ as mmtk_version
    mmtk_minversion = LooseVersion(vstring = '2.6.1')
    mmtk_version = LooseVersion(vstring = mmtk_version)
    
    if cmpver(mmtk_version, mmtk_minversion, pad = 0) == -1: raise 

except :
    print 'nMOLDYN needs MMTK %s or higher. Installation stopped.' % mmtk_minversion
    raise SystemExit

# Check that we have the minimal PyWin32 version.
if sys.platform == 'win32':
    try:
        pywin32_minversion = LooseVersion(vstring = '210')
        pywin32file = os.path.join(python_site_packages,'pywin32.version.txt')
        pywin32_version = str(file(pywin32file, 'r').read())
        
        if cmpver(pywin32_version, pywin32_minversion, pad = 0) == -1: raise 

    except :
        print 'nMOLDYN needs PyWin32 %s or higher. Installation stopped.' % pywin32_minversion
        raise SystemExit

# The nMOLDYN packages
packages = ['nMOLDYN',\
            'nMOLDYN.Analysis',\
            'nMOLDYN.Chemistry',\
            'nMOLDYN.Core',\
            'nMOLDYN.GUI',\
            'nMOLDYN.GUI.HTMLReader',\
            'nMOLDYN.Mathematics',\
            'nMOLDYN.Patches',\
            'nMOLDYN.Tests',\
            'nMOLDYN.Tests.ADOS',\
            'nMOLDYN.Tests.ARA',\
            'nMOLDYN.Tests.ARDCSF',\
            'nMOLDYN.Tests.ARDISF',\
            'nMOLDYN.Tests.AVACF',\
            'nMOLDYN.Tests.CDOS',\
            'nMOLDYN.Tests.CVACF',\
            'nMOLDYN.Tests.DCSF',\
            'nMOLDYN.Tests.DISF',\
            'nMOLDYN.Tests.DISFGA',\
            'nMOLDYN.Tests.EISF',\
            'nMOLDYN.Tests.MSD',\
            'nMOLDYN.Tests.RCF',\
            'nMOLDYN.Tests.nMOLDYN_Reference',\
            'nMOLDYN.Tests.SerialVsParallel',\
            ]

# The paths used to set the package datas.
paths = [
    os.path.join('Database'),\
    os.path.join('Database','Atoms'),\
    os.path.join('Database','Analysis'),\
    os.path.join('Doc'),\
    os.path.join('Doc','UsersGuide','HTML'),\
    os.path.join('Doc','API','HTML'),\
    os.path.join('Images'),\
    os.path.join('Patches'),\
    os.path.join('Trajectories'),\
]

if sys.platform == 'win32':
    paths.append(os.path.join('win32_files'))
    
# The package_data and other additional files.
# The additional data files.
options = {}
scripts = ['nMOLDYNStart.py']    
data_files = []
package_data = {'nMOLDYN': []}
for dir in paths:
    for f in glob(os.path.join('nMOLDYN', dir, '*')):
        if (f[-3:] != '.py') and (f[-4:-1] != '.py') and (os.path.isfile(f)):
            package_data['nMOLDYN'].append(os.path.join(dir,os.path.basename(f)))
package_data['nMOLDYN'].append(os.path.join('Doc','UsersGuide','Hyperlinks.xml'))
package_data['nMOLDYN'].append(os.path.join('Doc','UsersGuide','PDF','nMOLDYN_ug.pdf'))
package_data['nMOLDYN'].append(os.path.join('Doc','API','PDF','api.pdf'))

if sys.platform == 'win32':    
    data_files.append(('DLLs',[os.path.join('Lib', 'netcdf.dll')]))
    
    scripts.append('nMOLDYN_win32_postinstall.py')
    scripts.append('nMOLDYN_win32_uninstall.py')
    
    options['bdist_wininst'] = {'install_script' : 'nMOLDYN_win32_postinstall.py',\
                                'bitmap' : os.path.join('nMOLDYN', 'win32_files', 'nMOLDYN.bmp')}
    
requiredModules=['Python (>= 2.4)',\
                 'numpy (>= 1.2)',\
                 'matplotlib (>= 0.98)',\
                 'Pyro (>= 3.9)',\
                 'Scientific (>= 2.8)',\
                 'MMTK (>= 2.6.1)']

if sys.platform == 'win32':
    requiredModules.append('Pywin32 (>= 210.0)')	
    
# System-specific optimization options
low_opt = []
if sys.platform != 'win32' and sysconfig['CC'][:3] == 'gcc':
    low_opt = ['-O0']
    
high_opt = []
if sys.platform[:5] == 'linux' and sysconfig['CC'][:3] == 'gcc':
    high_opt = ['-O3', '-ffast-math', '-fomit-frame-pointer', '-fkeep-inline-functions']
    
if sys.platform == 'darwin' and sysconfig['CC'][:3] == 'gcc':
    high_opt = ['-O3', '-ffast-math', '-fomit-frame-pointer', '-fkeep-inline-functions', '-falign-loops=16']
high_opt.append('-g')

compile_args = []
compile_args.extend(low_opt)
compile_args.extend(high_opt)

link_args = []

define_macros = [("NUMPY", None)]
if platform.system() == "Windows" and platform.architecture()[0] == "64bit":
    import shutil
    dll = 'python%s.dll' % (''.join(platform.python_version_tuple()[:2]),)
    pyDLLPath = os.path.join(sys.prefix,"DLLs") 
    shutil.copyfile(os.path.join(os.environ["WINDIR"],"system32",dll), pyDLLPath)
    library_dirs = [pyDLLPath]
    define_macros.append(("MS_WIN64",None))
    link_args.append("-Wl,--allow-multiple-definition")    
else:
    library_dirs = []

include_dirs = [numpy.get_include(), os.path.abspath(os.path.join('Include','nMOLDYN'))]

headers = glob(os.path.join ('Include', 'nMOLDYN', '*.h'))

setup (name = "nMOLDYN",
       version = pkginfo.__version__,\
       description = "Analysis of Molecular Dynamics trajectories",\
       long_description =
"""nMOLDYN is an interactive program for the analysis of Molecular
Dynamics simulations. It is especially designed for the computation
and decomposition of neutron scattering spectra. The structure and
dynamics of the simulated systems can be characterized in terms of
various space and time correlation functions. To analyze the dynamics
of complex systems, rigid-body motions of arbitrarily chosen molecular
subunits can be studied.
""",\
       author = "E. Pellegrini, K. Hinsen, G.R. Kneller",\
       author_email = "pellegrini@ill.fr, research@khinsen.fastmail.net, kneller@cnrs-orleans.fr",\
       maintainer = "E. Pellegrini, K. Hinsen",\
       maintainer_email = "pellegrini@ill.fr, research@khinsen.fastmail.net",\
       url = "http://dirac.cnrs-orleans.fr/nMOLDYN/",\
       license = "CeCILL",\
       packages = packages,\
       package_data = package_data,\
       headers = headers,\
       platforms = ['Unix','Windows'],\
       ext_package = 'nMOLDYN.'+sys.platform,\
       ext_modules=[Extension("coordination_number",\
                              [os.path.join("Src","coordination_number.c")],\
                              define_macros = define_macros,\
                              library_dirs = library_dirs,\
                              extra_link_args = link_args,\
                              extra_compile_args = compile_args,\
                              include_dirs = include_dirs,\
                              ),\
                    Extension("distance_histogram",\
                              [os.path.join("Src","distance_histogram.c")],\
                              define_macros = define_macros,\
                              library_dirs = library_dirs,\
                              extra_link_args = link_args,\
                              extra_compile_args = compile_args,\
                              include_dirs = include_dirs,\
                              ),\
                    Extension("order_parameter",\
                              [os.path.join("Src","order_parameter.c")],\
                              define_macros = define_macros,\
                              library_dirs = library_dirs,\
                              extra_link_args = link_args,\
                              extra_compile_args = compile_args,\
                              include_dirs = include_dirs,\
                              ),\
                    Extension("smoothed_static_coherent_structure_factor",\
                              [os.path.join("Src","smoothed_static_coherent_structure_factor.c")],\
                              define_macros = define_macros,\
                              library_dirs = library_dirs,\
                              extra_link_args = link_args,\
                              extra_compile_args = compile_args,\
                              include_dirs = include_dirs,\
                              ),\
                    Extension("hbond_detection",\
                              [os.path.join("Src","hbond_detection.c")],\
                              define_macros = define_macros,\
                              library_dirs = library_dirs,\
                              extra_link_args = link_args,\
                              extra_compile_args = compile_args,\
                              include_dirs = include_dirs,\
                              )],\
       data_files = data_files,\
       requires = requiredModules,\
       scripts = scripts,\
       options = options
       )
