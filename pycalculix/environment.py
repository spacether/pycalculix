"""This module sets the dpi and the paths to gmsh, ccx, and cgx.

Attributes:
    DPI (None or float): if high dpi windows 8 monitor value is set, otherwise
        value is None.
    CCX (str): path to Calculix ccx, the sovler.
    CGX (str): path to Calculix cgx, the preprocessor/postprocessor/mesher.
    GMSH (str): path to the Gmsh mesher.
"""

import sys # needed to check if 32 or 64 bit interpreter
import os # used to delete files written by cgx
import platform # need this to check for win and do the below dpi fix
import ctypes # needed to see if it's 32-bit and to get correct win version

# http://stackoverflow.com/questions/19128219/detect-windows-8-1-in-python/22325767#22325767
class OSVERSIONINFOEXW(ctypes.Structure):
    """Returns object w/ attributes that will identify window os version"""
    _fields_ = [('dwOSVersionInfoSize', ctypes.c_ulong),
                ('dwMajorVersion', ctypes.c_ulong),
                ('dwMinorVersion', ctypes.c_ulong),
                ('dwBuildNumber', ctypes.c_ulong),
                ('dwPlatformId', ctypes.c_ulong),
                ('szCSDVersion', ctypes.c_wchar*128),
                ('wServicePackMajor', ctypes.c_ushort),
                ('wServicePackMinor', ctypes.c_ushort),
                ('wSuiteMask', ctypes.c_ushort),
                ('wProductType', ctypes.c_byte),
                ('wReserved', ctypes.c_byte)]

def get_version():
    """Get's the OS version. Returns a float of OS_MAJOR.OS_MINOR
    """
    os_version = OSVERSIONINFOEXW()
    os_version.dwOSVersionInfoSize = ctypes.sizeof(os_version)
    retcode = ctypes.windll.Ntdll.RtlGetVersion(ctypes.byref(os_version))
    if retcode != 0:
        print("Failed to get OS version")
        return -1.0

    major = str(os_version.dwMajorVersion)
    minor = str(os_version.dwMinorVersion)
    return float(major+'.'+minor)

def get_dpi():
    """Returns an int of current DPI for windows or None."""
    if platform.system() == 'Windows':
        version = get_version()
        if version >= 6.2:
            # Windows 8.0 is 6.2.9200
            # Only run on windows 8 or higher to check high dpi monitors

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
        #mysys = 'linux'+bits
        # assume that the user has the programs installed separately
        ccx = 'ccx'
        cgx = 'cgx'
        gmsh = 'gmsh'
    return [ccx, cgx, gmsh]

DPI = get_dpi()
CCX, CGX, GMSH = get_paths()
