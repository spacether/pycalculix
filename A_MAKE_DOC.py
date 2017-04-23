"""
This makes a docs folder which contains html documentation for the project.
"""
import shutil   # needed for folder removel(cleanum)
import subprocess   # needed to call python setup.py to make distribution
import os   #check if folder exists

distrib = 'pycalculix'
runbool = False

if __name__ == "__main__":
    runstr = input('Make the documentation? (Y/y) ')
    if runstr in ['Y','y']:
        runbool = True
else:
    runbool = True

if runbool:
    # remove folders before making them
    folders_to_rem = ['docs','../pycalculix_docs']  
    for fold in folders_to_rem:
        if os.path.isdir(fold):
            print ('Deleting folder: '+fold)
            shutil.rmtree(fold)

    # call sphinx-apidoc which auto-makes the scaffolding
    call_list = ['sphinx-apidoc', '-F',
                 '-H', 'pycalculix', '-A', 'Justin Black', '-V', '0.9.3',
                 '-e', '-o', 'docs', 'pycalculix']
    subprocess.call(call_list, shell=True)
    
    # edit the conf.py file to include the google formatting plugin
    fname = 'docs/conf.py'
    lines = []
    with open(fname, 'r') as f:
        lines = f.readlines()
    for (ind, line) in enumerate(lines):
        if "    'sphinx.ext.autodoc'," in line:
            lines[ind] =  lines[ind].rstrip()+"'sphinxcontrib.napoleon',\n"
            pass
        elif "#sys.path" in line:
            folder = '..'
            lines[ind] = "sys.path.insert(0, os.path.abspath('"+folder+"'))\n"
    with open(fname, 'w') as f:
        f.writelines(lines)
    print('conf.py edited to include code location and napoleon extension.')
    print('Google docstrings will now be read by sphinx.')

    # call sphinx to make the documentation in the docs folder
    subprocess.call("sphinx-build -b html docs ../pycalculix_docs", shell=True)

    # make needed nojekyl file
    jfile = open('../pycalculix_docs/.nojekyll', 'w').close()

    # remove folders before making them
    folders_to_rem = ['../pycalculix_docs/.doctrees']  
    for fold in folders_to_rem:
        if os.path.isdir(fold):
            print ('Deleting folder: '+fold)
            shutil.rmtree(fold)

else:
    print ('User chose not to make documentation.')

if __name__ == "__main__":
    input('Press enter to exit 0_MAKE_DOC.py ')
