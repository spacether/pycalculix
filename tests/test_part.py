import glob
import os
import subprocess
import unittest

import pycalculix as pyc

class TestPart(unittest.TestCase):
    """Test part module"""
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

    def test_part_sides_perfect(self):
        # a hole with multiple lines for each side sets them correctly
        proj_name = 'test_part_sides'

        model = pyc.FeaModel(proj_name)
        model.set_units('m') # this sets dist units to meters
        part = pyc.Part(model)

        divot_width = 1
        corner_width = 2

        part.goto(0, 0)
        left_points = [
            (corner_width, 0),
            (corner_width + divot_width, divot_width),
            (corner_width + divot_width*2, 0),
            (corner_width*2 + divot_width*2, 0),
        ]
        top_points = [
            (corner_width*2 + divot_width*2, corner_width),
            (corner_width*2 + divot_width, corner_width + divot_width),
            (corner_width*2 + divot_width*2, corner_width + divot_width*2),
            (corner_width*2 + divot_width*2, corner_width*2 + divot_width*2),
        ]
        right_points = [
            (corner_width + divot_width*2, corner_width*2 + divot_width*2),
            (corner_width + divot_width, corner_width*2 + divot_width),
            (corner_width, corner_width*2 + divot_width*2),
            (0, corner_width*2 + divot_width*2),
        ]
        bottom_points = [
            (0, corner_width + divot_width*2),
            (divot_width, corner_width + divot_width),
            (0, corner_width),
            (0, 0),
        ]
        all_points = left_points + top_points + right_points + bottom_points
        for x, y in all_points:
            part.draw_line_to(x, y)

        self.assertEqual([line.get_name() for line in part.left], ['L0', 'L3'])
        self.assertEqual([line.get_name() for line in part.top], ['L4', 'L7'])
        self.assertEqual([line.get_name() for line in part.right], ['L8', 'L11'])
        self.assertEqual([line.get_name() for line in part.bottom], ['L12', 'L15'])

    def test_part_sides_skewed(self):
        # top right point in range, bottom right point out of range
        # left side perfet (2 found)
        # top in range but skewed on 2nd line (2 found)
        # right in range but skewed on 1st line, 2nd out of range (1 found)
        # bottom 1st out of range, in range but skewed on 2nd line, 2 (1 found)
        initial_accuracy_val = pyc.geometry.ACC
        pyc.geometry.ACC = 0.2
        offset_in_range = 0.1
        offset_too_big = 0.4

        proj_name = 'test_part_sides_skewed'

        model = pyc.FeaModel(proj_name)
        model.set_units('m') # this sets dist units to meters
        part = pyc.Part(model)

        divot_width = 1
        corner_width = 2

        part.goto(0, 0)
        left_points = [
            (corner_width, 0),
            (corner_width + divot_width, divot_width),
            (corner_width + divot_width*2, 0),
            (corner_width*2 + divot_width*2, 0),
        ]
        top_points = [
            (corner_width*2 + divot_width*2, corner_width),
            (corner_width*2 + divot_width, corner_width + divot_width),
            (corner_width*2 + divot_width*2, corner_width + divot_width*2),
            (
                corner_width*2 + divot_width*2 - offset_in_range,
                corner_width*2 + divot_width*2 - offset_in_range
            ),
        ]
        right_points = [
            (corner_width + divot_width*2, corner_width*2 + divot_width*2),
            (corner_width + divot_width, corner_width*2 + divot_width),
            (corner_width, corner_width*2 + divot_width*2),
            (0 + offset_too_big, corner_width*2 + divot_width*2 - offset_too_big),
        ]
        bottom_points = [
            (0, corner_width + divot_width*2),
            (divot_width, corner_width + divot_width),
            (0, corner_width),
            (0, 0),
        ]
        all_points = left_points + top_points + right_points + bottom_points
        for x, y in all_points:
            part.draw_line_to(x, y)

        self.assertEqual([line.get_name() for line in part.left], ['L0', 'L3'])
        self.assertEqual([line.get_name() for line in part.top], ['L4', 'L7'])
        self.assertEqual([line.get_name() for line in part.right], ['L8'])
        self.assertEqual([line.get_name() for line in part.bottom], ['L15'])
        pyc.geometry.ACC = initial_accuracy_val

if __name__ == '__main__':
    unittest.main()
