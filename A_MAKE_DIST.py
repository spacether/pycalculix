"""
This file makes distributable packages in the zip file format
It deletes folders before making them, then does cleanup after it's done

Results:
    dist folder
        Has the pyalculix zip file one installs with pip
    examples folder
        Has the example py files only
    pycalculix
        Source code for the distribution, and binaries of gmsh and calculix
"""
import shutil   # needed for folder removel(cleanum) and making zip file
import subprocess   # needed to call python setup.py to make distribution
import os   #check if folder exists

distrib = 'pycalculix'
runstr = input('Make the distribution? (Y/y) ')

if runstr in ['Y','y']:

    # remove folders before making them
    folders_to_rem = ['dist', '%s.egg-info' % (distrib)]    
    for fold in folders_to_rem:
        if os.path.isdir(fold):
            print ('Deleting folder: '+fold)
            shutil.rmtree(fold)
    
    # call python
    subprocess.call("python setup.py sdist --formats=zip", shell=True)
    
    # delete egg-info folder
    folders_to_rem = ['%s.egg-info' % (distrib)]
    for fold in folders_to_rem:
        if os.path.isdir(fold):
            print ('Deleting folder: '+fold)
            shutil.rmtree(fold)
    
    # clean up the examples folder
    folder_to_clean = 'examples'
    print('Cleaning folder: '+folder_to_clean)
    for f in os.listdir(folder_to_clean):
        f = os.path.join(folder_to_clean, f)
        ext = os.path.splitext(f)[1]
        if ext not in ['.py', '.pdf', '.dxf']:
            os.remove(f)

    # zip the examples directory
    fromdir = 'examples'
    tofile = 'dist/examples'
    shutil.make_archive(tofile, 'zip', fromdir)

    # make documentation, zip it and move it into the dist folder
    import A_MAKE_DOC as docmake
    fromdir = '../pycalculix_docs'
    tofile = 'dist/documentation'
    shutil.make_archive(tofile, 'zip', fromdir)
    
    # overwrite current installation?
    over = input('Overwrite the current installation? (Y/y) ')
    if over in ['Y','y']:
        # check if the package is installed
        import pip
        installed_packs = pip.get_installed_distributions()
        installed_packs = [ i.key for i in installed_packs ] #get names
        if distrib in installed_packs:
            print('Uninstalling...')
            runstr = "pip uninstall -y %s" % (distrib)
            subprocess.call(runstr, shell=True)
        else:
            print('The library has not been installed yet')
        
        #install the library
        print('Installing...')
        zipfile = os.listdir('dist')
        for fname in zipfile:
            if distrib in fname:
                zipfile = os.path.join('dist',fname)
                subprocess.call("pip install %s" % (zipfile), shell=True)        

    # remove __pycache__ folder
    folders_to_rem = ['__pycache__']    
    for fold in folders_to_rem:
        if os.path.isdir(fold):
            print ('Deleting folder: '+fold)
            shutil.rmtree(fold)

else:
    print ('User chose not to make distribution')
    
input('Press enter to exit 0_MAKE_DIST.py ')
    

