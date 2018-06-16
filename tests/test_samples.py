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

    def example_tester(self, file_name, args=['-tri', '-nogui']):
        command_str = 'python examples/%s %s' % (file_name,
                                                  ' '.join(args))
        output = subprocess.check_output(
            command_str, stderr=subprocess.STDOUT, shell=True,
            timeout=60,
            universal_newlines=True)
        # this raises subprocess.CalledProcessError if it fails

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

    def test_hole_quarter(self, file_name='hole-in-plate-quarter.py'):
        self.example_tester(file_name)

    def test_hole_kt(self, file_name='hole-kt-study.py'):
        self.example_tester(file_name)

    def test_import_dxf1(self, file_name='import-dxf-1.py'):
        self.example_tester(file_name)

    def test_import_dxf2(self, file_name='import-dxf-2.py'):
        self.example_tester(file_name)

    def test_line_loop(self, file_name='line-loop.py'):
        self.example_tester(file_name)

    def test_multihole(self, file_name='multihole.py'):
        self.example_tester(file_name)

    def test_pinned_plate(self, file_name='pinned-plate.py'):
        self.example_tester(file_name)

    def test_pipe_crush_elastic(self,
                                file_name='pipe-crush-elastic.py'):
        self.example_tester(file_name)

    def test_rounded_square_ccw(self,
                                file_name='rounded-square-ccw.py'):
        self.example_tester(file_name)

    def test_rounded_square_cw(self,
                               file_name='rounded-square-cw.py'):
        self.example_tester(file_name)

    # note, once tests are written, make sure to add travisci file too
    # to ensure that it works on mac and linux

if __name__ == '__main__':
    unittest.main()
