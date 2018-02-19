import glob
import os
import subprocess
import unittest

import pycalculix as pyc

class TestExamples(unittest.TestCase):
    """Run all examples and make sure no exceptions are thrown"""
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

    def example_tester(self, file_name):
        ex = None
        try:
            command_str = 'python3 examples/%s -nogui' % file_name
            output = subprocess.check_output(
                command_str, stderr=subprocess.STDOUT, shell=True,
                timeout=120,
                universal_newlines=True)
        except subprocess.CalledProcessError as ex:
            print("Status : FAIL", ex.returncode, ex.output)
        else:
            print("Output: \n{}\n".format(output))
        self.assertIsNone(ex)

    def test_compr_rotor(self, file_name='compr-rotor.py'):
        self.example_tester(file_name)

    def test_dam_eplot(self, file_name='dam-eplot.py'):
        self.example_tester(file_name)

    def test_dam_times(self, file_name='dam-times.py'):
        self.example_tester(file_name)

    def test_dam(self, file_name='dam.py'):
        self.example_tester(file_name)

    def test_hole_fancy(self, file_name='hole-fancy.py'):
        self.example_tester(file_name)

    def test_hole_full(self, file_name='hole-in-plate-full.py'):
        self.example_tester(file_name)

    def test_multihole(self, file_name='multihole.py'):
        self.example_tester(file_name)

if __name__ == '__main__':
    unittest.main()
