import glob
import os
import subprocess
import unittest

import pycalculix as pyc

class TestPart(unittest.TestCase):
    """Test meshing module"""
    def tearDown(self):
        extenstion_to_del = ['fbd',
                             'inp',
                             'geo',
                             'msh',
                             'frd',
                             'dat',
                             'png',
                             'cvg',
                             'sta',
                             'out']
        local_files = glob.glob('*.*') # files that have extensions
        for local_file in local_files:
            extension = local_file.split('.')[-1]
            if extension in extenstion_to_del:
                os.unlink(local_file)

    def test_part_sides(self):
        # a hole with multiple lines for each side sets them correctly
        pass

    def test_part_sides(self):
        # a hole with multiple lines for each side sets them correctly
        # even if one of them is skewed
        pass

if __name__ == '__main__':
    unittest.main()
