# Pour compiler:
#      python setup.py build_ext --inplace

from distutils.core import setup, Extension
from Cython.Distutils import build_ext

import os
import sys
import platform

import numpy

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

include_dirs = [numpy.get_include(), os.path.abspath(os.path.join('..', 'Include','nMOLDYN'))]

setup (name = "nMOLDYNPyrexModules",
       version = "0.1",

       ext_modules = [Extension('distance_histogram', ['distance_histogram.pyx'],\
                                define_macros = define_macros,\
                                library_dirs = library_dirs,\
                                extra_link_args = link_args,\
                                include_dirs = include_dirs,\
                                ),\
                      Extension('smoothed_static_coherent_structure_factor', ['smoothed_static_coherent_structure_factor.pyx'],\
                                define_macros = define_macros,\
                                library_dirs = library_dirs,\
                                extra_link_args = link_args,\
                                include_dirs = include_dirs,\
                                ),\
                      Extension('coordination_number', ['coordination_number.pyx'],\
                                define_macros = define_macros,\
                                library_dirs = library_dirs,\
                                extra_link_args = link_args,\
                                include_dirs = include_dirs,\
                                ),\
                      Extension('hbond_detection', ['hbond_detection.pyx'],\
                                define_macros = define_macros,\
                                library_dirs = library_dirs,\
                                extra_link_args = link_args,\
                                include_dirs = include_dirs,\
                                ),\
                      Extension('order_parameter', ['order_parameter.pyx'],\
                                define_macros = define_macros,\
                                library_dirs = library_dirs,\
                                extra_link_args = link_args,\
                                include_dirs = include_dirs,\
                                )],\
       cmdclass = {'build_ext': build_ext}
       )
