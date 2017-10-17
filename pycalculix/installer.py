import shutil
import shlex
import subprocess
import sys

from sys import platform as _platform

def add():
    osname = None
    is_64bit = sys.maxsize > 2**32
    bitsize_dict = {True: '64', False: '32'}
    bitsize = bitsize_dict[is_64bit]
    printos = lambda os, bits=bitsize: print('Detected {} {} bit'.format(os, bits))
    if _platform == "linux" or _platform == "linux2":
        printos('Linux')
    elif _platform == "darwin":
        printos('Mac OS X')
        brew_installed = shutil.which('brew')
        if not brew_installed:
            print('Installing brew')
            command_line = "/usr/bin/ruby -e \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)\""
            subprocess.check_call(command_line, shell=True)
        else:
            print('brew present')
        gmsh_installed = shutil.which('gmsh')
        if not gmsh_installed:
            print('Installing gmsh')
            command_line = "brew install gmsh"
            subprocess.check_call(command_line, shell=True)
        else:
            print('gmsh present')
        # need to fix ccx here
        ccx_installed = shutil.which('ccx')
        if not ccx_installed:
            print('Installing calculix (ccx)')
            command_line = "brew install homebrew/science/calculix-ccx"
            subprocess.check_call(command_line, shell=True)
            command = "ln -s /usr/local/bin/ccx_2.12 /usr/local/bin/ccx"
            subprocess.check_call(command_line, shell=True)
        else:
            print('calculix (ccx) present')
    elif _platform == "win32":
        printos('Windows')
    elif _platform == "win64":
        printos('Windows')

def remove():
    osname = None
    is_64bit = sys.maxsize > 2**32
    bitsize_dict = {True: '64', False: '32'}
    bitsize = bitsize_dict[is_64bit]
    printos = lambda os, bits=bitsize: print('Detected {} {} bit'.format(os, bits))
    if _platform == "linux" or _platform == "linux2":
        printos('Linux')
    elif _platform == "darwin":
        printos('Mac OS X')
        brew_installed = shutil.which('brew')
        if not brew_installed:
            print('Installing brew')
            command_line = "/usr/bin/ruby -e \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)\""
            subprocess.check_call(command_line, shell=True)
        else:
            print('brew present')
        ccx_installed = shutil.which('ccx')
        if not ccx_installed:
            print('calculix (ccx) is not on your system')
        else:
            print('Removing calculix (ccx)')
            command_line = "rm /usr/local/bin/ccx"
            subprocess.check_call(command_line, shell=True)
            command_line = "brew uninstall homebrew/science/calculix-ccx"
            subprocess.check_call(command_line, shell=True)
        gmsh_installed = shutil.which('gmsh')
        if not gmsh_installed:
            print('gmsh is not on your system')
        else:
            print('Removing gmsh')
            command_line = "brew uninstall gmsh"
            subprocess.check_call(command_line, shell=True)
        # need to fix ccx here
    elif _platform == "win32":
        printos('Windows')
    elif _platform == "win64":
        printos('Windows')
