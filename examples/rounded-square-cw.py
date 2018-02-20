#!/usr/bin/env python3
import sys

import pycalculix as pyc

model_name = 'rounded-square-cw'
model = pyc.FeaModel(model_name)
model.set_units('m')

# the below boolean sets whether or not to show gui plots
# testing passes in -nogui
show_gui = True
if len(sys.argv) == 2 and sys.argv[-1] == '-nogui':
    show_gui = False

length = 2
thickness = 1
radius = 0.5
radius_inner = 1
axial = 10

part = pyc.Part(model)
part.goto(radius_inner, axial)
part.draw_line_rad(thickness*0.6)
part.draw_line_to(radius_inner+thickness, axial)
line_1 = part.draw_line_ax(length)[0]
line_2 = part.draw_line_rad(-thickness)[0]
part.draw_line_ax(-length*0.5)
part.draw_line_to(radius_inner, axial)
part.fillet_lines(line_1, line_2, radius)

model.plot_geometry(model_name + '_geom', display=show_gui)
part.chunk()
model.plot_geometry(model_name + '_chunked', display=show_gui)
model.view.print_summary()

model.set_etype('axisym', part)
model.mesh(1.0, 'gmsh')
model.plot_elements(model_name+'_elements', display=show_gui)
model.view.print_summary()
model.print_summary()
