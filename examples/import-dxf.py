#!/usr/bin/env python3
import pycalculix as pyc

model_name = 'import-dxf'
model = pyc.FeaModel(model_name)
model.set_units('m')

#fname = 'test.dxf'
fname = 'kontrola.dxf'
importer = pyc.CadImporter(model, fname, swapxy=True)
parts = importer.load()
model.plot_geometry(model_name+'_imported')
#parts[0].chunk()
model.plot_geometry(model_name+'_chunked_areas', pnum=False, lnum=False)
model.plot_geometry(model_name+'_chunked_lines', anum=False, pnum=False)
model.plot_geometry(model_name+'_chunked_points', anum=False, lnum=False)

model.view.print_summary()

model.set_etype('axisym', parts)
model.set_eshape('quad', 2)
model.mesh(1.0, 'gmsh')
model.plot_elements(model_name+'_elements')
model.view.print_summary()
