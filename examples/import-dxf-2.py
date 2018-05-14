#!/usr/bin/env python3
import os
import sys

import pycalculix as pyc

model_name = 'import-dxf-2'
model = pyc.FeaModel(model_name)
model.set_units('m')

# set whether or not to show gui plots
show_gui = True
if '-nogui' in sys.argv:
    show_gui = False
# set element shape
eshape = 'quad'
if '-tri' in sys.argv:
    eshape = 'tri'

#fname = 'test.dxf'
abs_path = os.path.dirname(os.path.abspath(__file__))
fname = os.path.join(abs_path, '%s.dxf' % model_name)
importer = pyc.CadImporter(model, fname, swapxy=True)
parts = importer.load()
model.plot_geometry(model_name+'_imported', display=show_gui)
#parts[0].chunk()
model.plot_geometry(model_name+'_areas', pnum=False,
                    lnum=False, display=show_gui)
model.plot_geometry(model_name+'_lines', anum=False,
                    pnum=False, display=show_gui)
model.plot_geometry(model_name+'_points', anum=False,
                    lnum=False, display=show_gui)


model.view.print_summary()

model.set_etype('axisym', parts)
model.set_eshape(eshape, 2)
model.mesh(1.0, 'gmsh')
model.plot_elements(model_name+'_elements', display=show_gui)
model.view.print_summary()
