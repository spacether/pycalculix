"""
This file ceans the examples folder, leaving only .py files in it.
"""
import shutil   # needed for folder removel(cleanum)
import subprocess   # needed to call python setup.py to make distribution
import os   #check if folder exists

distrib = 'pycalculix'
runstr = input('Clean examples folder? (Y/y) ')

if runstr in ['Y','y']:

    # clean up the examples folder
    folder_to_clean = 'examples'
    print('Cleaning folder: '+folder_to_clean)
    for f in os.listdir(folder_to_clean):
        f = os.path.join(folder_to_clean, f)
        ext = os.path.splitext(f)[1]
        if ext not in ['.py', '.pdf', '.dxf']:
            os.remove(f)

else:
    print ('User chose not to clean the examples folder')
    
input('Press enter to exit 0_CLEAN_EXAMPLES.py ')
    

