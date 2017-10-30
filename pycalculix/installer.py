
import glob
import os
import re
import shutil
import shlex
import subprocess
import sys
import urllib
import zipfile

from sys import platform as platform

def add():
    osname = None
    is_64bit = sys.maxsize > 2**32
    bitsize_dict = {True: 64, False: 32}
    bitsize = bitsize_dict[is_64bit]
    printos = lambda os, bits=bitsize: print('Detected {} {} bit'.format(os, bits))
    if platform == "linux" or platform == "linux2":
        printos('Linux')
    elif platform == "darwin":
        printos('Mac OS X')
        mac_add()
    elif platform == "win32":
        printos('Windows')
        windows_add(bitsize)
    elif platform == "win64":
        printos('Windows')
        windows_add(bitsize)
    print('Done!')

def remove():
    osname = None
    is_64bit = sys.maxsize > 2**32
    bitsize_dict = {True: 64, False: 32}
    bitsize = bitsize_dict[is_64bit]
    printos = lambda os, bits=bitsize: print('Detected {} {} bit'.format(os, bits))
    if platform == "linux" or platform == "linux2":
        printos('Linux')
    elif platform == "darwin":
        printos('Mac OS X')
        mac_remove()
    elif platform == "win32":
        printos('Windows')
        windows_remove(bitsize)
    elif platform == "win64":
        printos('Windows')
        windows_remove(bitsize)
    print('Done!')

def mac_add():
    """Adds programs on mac, uses brew"""
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
    ccx_installed = shutil.which('ccx')
    if not ccx_installed:
        print('Installing calculix (ccx)')
        command_line = "brew install homebrew/science/calculix-ccx"
        subprocess.check_call(command_line, shell=True)
        # add search here to find the path to the original gmsh program
        command = "ln -s /usr/local/bin/ccx_2.12 /usr/local/bin/ccx"
        subprocess.check_call(command_line, shell=True)
    else:
        print('calculix (ccx) present')

def mac_remove():
    """Removes programs on mac, uses brew"""
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

def windows_add(bitsize):
    """Adds programs on windows"""
    gmsh_installed = shutil.which('gmsh')
    if not gmsh_installed:
        print('Installing gmsh')
        win_add_from_url(bitsize, 'http://gmsh.info/bin/Windows/', 'gmsh')
    else:
        print('gmsh present')

    ccx_installed = shutil.which('ccx')
    if not ccx_installed:
        print('Installing calculix (ccx)')
        win_add_ccx(bitsize, "https://sourceforge.net/projects/calculixforwin/files/03.2/", "ccx")
    else:
        print('calculix (ccx) present')

    # ccx Notes
    # to get 32 bit and 64 bit can use either:
    # http://www.bconverged.com/data/content/CalculiX_2_7_win_003.zip
    # with this, one has to manually pick either 32 or 64bit
    # eliminated
    # then I need to add the folder: C:\Program Files (x86)\bConverged\CalculiX\ccx to path

    # https://downloads.sourceforge.net/project/calculixforwin/03.2/CL32-win32bit.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fcalculixforwin%2Ffiles%2F03.2%2F&ts=1509080746&use_mirror=astuteinternet
    # from that unzipped file I can move ccx out from C:\Users\Justin\Downloads\CL32-win64bit\CL32-win64\bin\ccx
    # but keep the dlls too
    # use mklink or os.symlink to create a link to ccx

    # sourceforge also has a 64 bit version
    # github does not have a 32 bit version of the binary sadly

def windows_remove(bitsize):
    """Removes programs on windows"""
    gmsh_installed = shutil.which('gmsh')
    if not gmsh_installed:
        print('gmsh is not on your system')
    else:
        print('Removing gmsh')
        env_path = os.getenv('VIRTUAL_ENV', sys.exec_prefix)
        scripts_folder = '%s\Scripts\\' % env_path
        remove_like(scripts_folder, 'gmsh')

def remove_like(search_path, name):
    """
    Searches for files or folders matching name in search_path and deletes them
    """
    match_str = '%s*%s*' % (search_path, name)
    paths = glob.glob(match_str)
    for path in paths:
        print('Removing %s' % path)
        if os.path.isfile(path):
            os.unlink(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)

def win_add_from_url(bitsize, binaries_url, program_name):
    # Needs the user agent header to exist for it to send back the content
    user_agent = ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:56.0) '
                  'Gecko/20100101 Firefox/56.0')
    headers = {'User-Agent': user_agent}

    zipfile_regex = '.*%s.*' % program_name
    url_choices = zipfile_by_bitsize(binaries_url, headers, zipfile_regex)

    zipfile_name = url_choices[bitsize]
    zipfile_folder_name = zipfile_name.split('.')[0]
    zipfile_url = binaries_url + zipfile_name
    print('Downloading %s from %s' % (program_name, zipfile_url))
    req = urllib.request.Request(zipfile_url, None, headers)
    with urllib.request.urlopen(req) as in_stream, open(zipfile_name, 'wb') as out_file:
        shutil.copyfileobj(in_stream, out_file)
    print('Unzipping %s' % program_name)
    zip_ref = zipfile.ZipFile(zipfile_name, 'r')
    zip_ref.extractall(None)
    zip_ref.close()
    print('Removing %s zipfile' % program_name)
    os.remove(zipfile_name)

    env_path = os.getenv('VIRTUAL_ENV', sys.exec_prefix)
    folder_to = '%s\Scripts\%s' % (env_path, zipfile_folder_name)
    exe_loc = '%s\%s.exe' % (folder_to, program_name)
    exe_link = '%s\Scripts\%s.exe' % (env_path, program_name)
    paths = [exe_link, folder_to]
    for path in  paths:
        if os.path.isfile(path):
            os.unlink(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)
    command_line = "move %s %s" % (zipfile_folder_name, folder_to)
    print('Installing %s to %s' % (program_name, folder_to))
    subprocess.check_call(command_line, shell=True)
    os.link(exe_loc, exe_link)

def zipfile_by_bitsize(binaries_url, headers, zipfile_regex):
    req = urllib.request.Request(binaries_url, None, headers)
    html = urllib.request.urlopen(req).read().decode('utf-8')
    urls = re.findall(r'href=[\'"]?([^\'" >]+)', html)
    pattern = re.compile(zipfile_regex)
    urls = [url for url in urls if pattern.match(url)]
    urls = urls[-2:]
    url_choices = {32: urls[0], 64: urls[1]}
    if 'win32' in urls[1] or 'Windows32' in urls[1]:
        url_choices = {32: urls[1], 64: urls[0]}
    return url_choices

def get_direct_url(source_page, headers):
    req = urllib.request.Request(source_page, None, headers)
    html = urllib.request.urlopen(req).read().decode('utf-8', 'ignore')
    html = html.replace("\"", "'")
    link_text_pos = html.find('direct link')
    href_pos = html[:link_text_pos].rfind('href')
    first_char = html.find("'", href_pos)+1
    last_quote = html.find("'", first_char)
    return html[first_char:last_quote]

def win_add_ccx(bitsize, binaries_url, program_name):
    # Needs the user agent header to exist for it to send back the content
    user_agent = ('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:56.0) '
                  'Gecko/20100101 Firefox/56.0')
    headers = {'User-Agent': user_agent}

    zipfile_regex = '.+CL.+win.+zip/download$'
    url_choices = zipfile_by_bitsize(binaries_url, headers, zipfile_regex)
    zipfile_name = url_choices[bitsize].split('/')[-2]

    zipfile_webpage_url = url_choices[bitsize]
    print(zipfile_webpage_url)
    zipfile_url = get_direct_url(zipfile_webpage_url, headers)
    print(zipfile_url)

    print('Downloading %s from %s' % (zipfile_name, zipfile_url))
    req = urllib.request.Request(zipfile_url, None, headers)
    with urllib.request.urlopen(req) as in_stream, open(zipfile_name, 'wb') as out_file:
        shutil.copyfileobj(in_stream, out_file)

    print('Unzipping %s' % program_name)
    zip_ref = zipfile.ZipFile(zipfile_name, 'r')
    zip_ref.extractall(None)
    zip_ref.close()
    print('Removing %s zipfile' % program_name)
    os.remove(zipfile_name)
