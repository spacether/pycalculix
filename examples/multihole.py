#!/usr/bin/env python3
import pycalculix as pyc

model_name = 'multihole'
model = pyc.FeaModel(model_name)
model.set_units('m')

width = 5
length = 8
radius = 0.5


part = pyc.Part(model)
part.goto(0, 0)
part.draw_line_rad(length)
part.draw_line_ax(width)
part.draw_line_rad(-length)
part.draw_line_ax(-width)

"""
# Circular Holes
segments = 4
radius = 0.5
hole_lines = part.draw_hole(2, 2, radius, segments)
model.set_ediv(hole_lines, 10)
hole_lines = part.draw_hole(4, 2, radius, segments)
model.set_ediv(hole_lines, 10)
hole_lines = part.draw_hole(6, 2, radius, segments)
model.set_ediv(hole_lines, 10)
"""

# Square Holes
hole_width = 1
hole_length = 1
part.goto(1, 1, holemode=True)
part.draw_line_rad(hole_length)
part.draw_line_ax(hole_width)
part.draw_line_rad(-hole_length)
part.draw_line_ax(-hole_width)

part.goto(1.5, 3, holemode=True)
part.draw_line_rad(hole_length)
part.draw_line_ax(hole_width)
part.draw_line_rad(-hole_length)
part.draw_line_ax(-hole_width)


model.plot_geometry(model_name+'_prechunk_areas', lnum=False, pnum=False)
model.plot_geometry(model_name+'_prechunk_lines', pnum=False)
model.plot_geometry(model_name+'_prechunk_points', lnum=False)

part.chunk(debug=[0,0])
model.plot_geometry(model_name+'_chunked_areas', lnum=False, pnum=False)

model.set_eshape('quad', 2)
model.set_etype('plstrain', part, 0.1)
model.mesh(0.7, 'gmsh')
model.plot_elements()

model.view.print_summary()
