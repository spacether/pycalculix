"""This module sets the dpi and the paths to gmsh, ccx, and cgx."""

import sys # needed to check if 32 or 64 bit interpreter
import os # used to delete files written by cgx
import platform # need this to check for win and do the below dpi fix

def get_dpi():
    """Returns an int of current DPI for windows or None."""
    if platform.system() == 'Windows':
        if float(platform.release()) >= 8:
            # Only run on windows 8 or higher to check high dpi monitors

            import ctypes
            user32 = ctypes.windll.user32
            w_curr = user32.GetSystemMetrics(0)
            user32.SetProcessDPIAware()
            w_phys = user32.GetSystemMetrics(0)
            curr_dpi = round(w_phys*96/w_curr, 0)

            from pylab import rcParams
            rcParams['figure.dpi'] = curr_dpi
            return curr_dpi
    else:
        return None

def get_paths():
    """Returns a list of paths to: [ccx, cgx, gmsh]."""
    # this points the program to the installed location of ccx, cgx, and gmsh
    bits = '32'
    mydir = os.path.dirname(__file__)
    (ccx, cgx, gmsh) = ('', '', '')
    if sys.maxsize > 2**32:
        # we are using a 64 bit interpreter
        bits = '64'
    if platform.system() == 'Windows':
        mysys = 'win'+bits
        ccx = os.path.join(mydir, 'calculix_'+mysys, 'ccx.exe')
        cgx = os.path.join(mydir, 'calculix_'+mysys, 'cgx.exe')
        gmsh = os.path.join(mydir, 'gmsh_win32', 'gmsh.exe')
    elif platform.system() == 'Linux':
        mysys = 'linux'+bits
        ccx = 'ccx'
        cgx = os.path.join(mydir, 'calculix_'+mysys, 'cgx')
        gmsh = 'gmsh'
    return [ccx, cgx, gmsh]

DPI = get_dpi()
[CCX, CGX, GMSH] = get_paths()
