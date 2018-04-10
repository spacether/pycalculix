import glob
import os
import subprocess
import unittest

import pycalculix as pyc

class TestMeshing(unittest.TestCase):
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

    def most_are_type(self, elements, el_shape):
        min_num_nodes = 3
        if el_shape == 'quad':
            min_num_nodes = 4
        # require 90% be asked type
        passing_qty = 0.9*len(elements)
        type_count = 0
        for element in elements:
            if len(element.nodes) % min_num_nodes == 0:
                type_count += 1
        print('Percent %ss %f' %
              (el_shape, 100*type_count/len(elements)))
        if type_count >= passing_qty:
            return True
        return False

    def test_tri_mesh(self, el_shape='tri'):
        proj_name = 'sample'
        length = 4.0 # hole diam

        model = pyc.FeaModel(proj_name)
        model.set_units('m') # this sets dist units to meters
        part = pyc.Part(model)

        part.goto(0, 0)
        part.draw_line_ax(length)
        part.draw_line_rad(length)
        part.draw_line_ax(-length)
        part.draw_line_rad(-length)

        # set the element type and mesh database
        model.set_eshape(el_shape, 2)
        model.mesh(1.0, 'gmsh') # mesh 1.0 fineness, smaller is finer
        model.plot_elements(proj_name+'_elem')   # plot part elements
        self.assertTrue(self.most_are_type(part.elements, el_shape))

    def test_quad_mesh(self, el_shape='quad'):
        proj_name = 'sample'
        length = 4.0 # hole diam

        model = pyc.FeaModel(proj_name)
        model.set_units('m') # this sets dist units to meters
        part = pyc.Part(model)

        part.goto(0, 0)
        part.draw_line_ax(length)
        part.draw_line_rad(length)
        part.draw_line_ax(-length)
        part.draw_line_rad(-length)

        # set the element type and mesh database
        model.set_eshape(el_shape, 2)
        model.mesh(1.0, 'gmsh') # mesh 1.0 fineness, smaller is finer
        model.plot_elements(proj_name+'_elem')   # plot part elements
        self.assertTrue(self.most_are_type(part.elements, el_shape))

    # add tests of fineness here

if __name__ == '__main__':
    unittest.main()
