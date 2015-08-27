# Add Python DLL directory to PATH

from distutils.version import LooseVersion
import os
import shutil
import sys
import time
import nMOLDYN    

def apply_patch_for_scientific_module():
    """Apply a patch for one file that is buggy in case of a win32 installation.
    """
    
    scientific_minversion = '2.8.0'
    import Scientific
    scientific_version = LooseVersion(vstring = Scientific.__version__)
    
    # Versions > 2.8.0 should not need the patch.
    if scientific_version > LooseVersion(vstring = scientific_minversion): return
    
    # The path to the ScientificPythin installation.
    scientific_path = Scientific.__path__[0]
    
    # The path to its DistributedComputing path.
    dir = os.path.join(scientific_path, 'DistributedComputing')
    
    # The file that will have to be replaced.
    src_file = os.path.join(dir, 'MasterSlave.py')
    
    # It will be copy first under this name.
    dst_file = os.path.join(dir, 'MasterSlave_copy_%s.py' % time.strftime("%a_%d_%b_%Y", time.localtime()))
    
    # The copy is done.
    shutil.copyfile(src_file, dst_file)
    
    # Some info.
    print 'The file %s was copied to %s before applying the patch.' % (src_file, dst_file)
        
    # The path to the patch.
    patch_file = os.path.join(nMOLDYN.__path__[0],'Patches','MasterSlave_patch.py')

    # Try to aply the patch by replacing the MasterSlave.py file by the patch.
    try:    
        shutil.copyfile(patch_file, src_file)
    # If this can not be done, this might be due to a permission problem.
    except:
        print 'The patch %s could not be applied. Perhaps a permission problem.'
        
    # If it is OK, print some info.
    else:    
        print 'The patch %s was successfully applied.' % patch_file

def install():
    """Will performes post-installation opeartions.
    """

    # First apply the ScientifcPython patch to make PyRo work on Windows.
    apply_patch_for_scientific_module()

    # The path to the nMOLDYN icon.
    launch_ico_path = os.path.join(nMOLDYN.__path__[0],'win32_files','nMOLDYN_Logo.ico')

    # The path to the nMOLDYN uninstall icon.
    uninstall_ico_path = os.path.join(nMOLDYN.__path__[0],'win32_files','uninstall_nMOLDYN.ico')
    
    # The path to the nMOLDYN starting script.
    script_path = os.path.join(sys.prefix, 'Scripts', 'nMOLDYNStart.py')
    
    # The path to the nMOLDYN uninstall script.
    uninstall_script_path = os.path.join(sys.prefix, 'Scripts', 'nMOLDYN_win32_uninstall.py')

    try:
        desktop_path = get_special_folder_path("CSIDL_COMMON_DESKTOPDIRECTORY")
        
    except OSError:
        desktop_path = get_special_folder_path("CSIDL_DESKTOPDIRECTORY")

    print 'Desktop path created in %s' % desktop_path
        
    desktop_shortcut = os.path.join(desktop_path, 'nMOLDYN 3.lnk')

    # Creates the desktop shortcut.                
    create_shortcut(script_path,\
                    'This is the shortcut for the nMOLDYN 3 launcher',\
                    desktop_shortcut,\
                    '',\
                    '',\
                    launch_ico_path,\
                    0)

    # Register that path with the uninstaller, so that it will be removed when the distribution is uninstalled.
    file_created(desktop_shortcut)
    
    try:
        start_path = get_special_folder_path("CSIDL_COMMON_PROGRAMS")
        
    except OSError:
        start_path = get_special_folder_path("CSIDL_PROGRAMS")
        
    print 'Start path created in %s' % start_path
            
    programs_path = os.path.join(start_path, "nMOLDYN 3")

    if not os.path.exists(programs_path):
        try :
            os.mkdir(programs_path)

        except OSError:
            pass

    # Register that directory with the uninstaller, so that it will be removed when the distribution is uninstalled.
    directory_created(programs_path)
    
    start_shortcut = os.path.join(programs_path, 'nMOLDYN 3.lnk')
    
    # Creates the Start Menu launching shortcut
    create_shortcut(script_path,\
                    'nMOLDYN 3 launcher',\
                    start_shortcut,\
                    '',\
                    '',\
                    launch_ico_path,\
                    0)
                    
    # Register that path with the uninstaller, so that it will be removed when the distribution is uninstalled.
    file_created(start_shortcut)

    start_uninstall_shortcut = os.path.join(programs_path, 'uninstall nMOLDYN 3.lnk')

    # Creates the Start Menu uninstall shortcut.
    create_shortcut(uninstall_script_path,\
                    'Uninstall nMOLDYN 3',\
                    start_uninstall_shortcut,\
                    '',\
                    '',\
                    uninstall_ico_path,\
                    0)
                    
    # Register that path with the uninstaller, so that it will be removed when the distribution is uninstalled.
    file_created(start_uninstall_shortcut)

if __name__=='__main__':
    if len(sys.argv) == 2 and sys.argv[1] == '-install':
        install()

