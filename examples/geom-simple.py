#!/usr/bin/env python3
import pycalculix as pyc

# Vertical hole in plate model, make model
proj_name = 'geom-simple'
model = pyc.FeaModel(proj_name)
model.set_units('m') # this sets dist units to meters

# Define variables we'll use to draw part geometry
diam = 2.000 # hole diam
width = 4.444 # plate width
ratio = diam/width
print('D=%f, H=%f, D/H=%f' % (diam, width, ratio))
length = 2*width  #plate length
width_quarter = width/2
length_quarter = length/2
rad = diam*0.5

# Draw part geometry, x, y = radial, axial
part = pyc.Part(model)
part.goto(0.0, rad)
part.draw_arc_angle(90, 0.0, 0.0)
part.draw_line_to(length_quarter, 0.0)
part.draw_line_ax(width_quarter)
part.draw_line_rad(-length_quarter*0.5)
part.draw_line_rad(-length_quarter*0.5)
part.draw_line_to(0.0, rad)
model.plot_lines(proj_name+'_lines') # view the geometry
