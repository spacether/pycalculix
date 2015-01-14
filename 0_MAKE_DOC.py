"""
This makes a docs folder which contains html documentation for the project.
"""
import shutil   # needed for folder removel(cleanum)
import subprocess   # needed to call python setup.py to make distribution
import os   #check if folder exists

distrib = 'pycalculix'
runstr = input('Make the documentation? (Y/y) ')

if runstr in ['Y','y']:

    # remove folders before making them
    folders_to_rem = ['docs_build','docs']  
    for fold in folders_to_rem:
        if os.path.isdir(fold):
            print ('Deleting folder: '+fold)
            shutil.rmtree(fold)

    # call sphinx-apidoc which auto-makes the scaffolding
    subprocess.call("sphinx-apidoc -F -e -o docs_build pycalculix", shell=True)
    
    # edit the conf.py file to include the google formatting plugin
    fname = 'docs_build/conf.py'
    lines = []
    with open(fname, 'r') as f:
        lines = f.readlines()
    for (ind, line) in enumerate(lines):
        if "    'sphinx.ext.autodoc'," in line:
            lines[ind] =  lines[ind].rstrip()+"'sphinxcontrib.napoleon',\n"
    with open(fname, 'w') as f:
        f.writelines(lines)
    print('conf.py edited to include the napoleon extension.')
    print('Google docstrings will now be read by sphinx.')

    # call sphinx to make the documentation in the docs folder
    subprocess.call("sphinx-build -b html docs_build docs", shell=True)

    # remove the docs_build folder
    folders_to_rem = ['docs_build']  
    for fold in folders_to_rem:
        if os.path.isdir(fold):
            print ('Deleting folder: '+fold)
            shutil.rmtree(fold)

else:
    print ('User chose not to make documentation.')
    
input('Press enter to exit 0_MAKE_DOC.py ')
    

