#!/usr/bin/env python3
import os
import sys

import pycalculix as pyc

model_name = 'import-dxf'
model = pyc.FeaModel(model_name)
model.set_units('m')

# the below boolean sets whether or not to show gui plots
# testing passes in -nogui
show_gui = True
if len(sys.argv) == 2 and sys.argv[-1] == '-nogui':
    show_gui = False

#fname = 'test.dxf'
abs_path = os.path.dirname(os.path.abspath(__file__))
fname = os.path.join(abs_path, 'kontrola.dxf')
importer = pyc.CadImporter(model, fname, swapxy=True)
parts = importer.load()
model.plot_geometry(model_name+'_imported')
#parts[0].chunk()
model.plot_geometry(model_name+'_chunked_areas', pnum=False,
                    lnum=False, display=show_gui)
model.plot_geometry(model_name+'_chunked_lines', anum=False,
                    pnum=False, display=show_gui)
model.plot_geometry(model_name+'_chunked_points', anum=False,
                    lnum=False, display=show_gui)

model.view.print_summary()

model.set_etype('axisym', parts)
model.set_eshape('quad', 2)
model.mesh(1.0, 'gmsh')
model.plot_elements(model_name+'_elements', display=show_gui)
model.view.print_summary()
