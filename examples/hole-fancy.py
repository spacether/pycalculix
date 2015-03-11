import pycalculix as pyc

# arc test sample
proj_name = 'hole-fancy'
model = pyc.FeaModel(proj_name)
model.set_units('m')    # this sets dist units to meters, labels our consistent units

# Define variables we'll use to draw part geometry
length = 8
width = 6
radius = 1
hole_width = width - 2*radius
hole_length = length - 2*radius

# Draw part geometry, you must draw the part CLOCKWISE, x, y = radial, axial
part = pyc.Part(model)
part.goto(length*0.5, -width*0.5)
part.draw_line_ax(width)
part.draw_line_rad(-length)
part.draw_line_ax(-width)
part.draw_line_rad(length)

# make hole
part.goto(hole_length*0.5, hole_width*0.5, holemode=True)
l_top = part.draw_line_ax(-hole_width)[0]
l_left = part.draw_line_rad(-hole_length)[0]
l_bot = part.draw_line_ax(hole_width)[0]
l_right = part.draw_line_rad(hole_length)[0]
arcs = []
arcs.append(part.fillet_lines(l_top, l_left, radius)[0])
arcs.append(part.fillet_lines(l_left, l_bot, radius)[0])
arcs.append(part.fillet_lines(l_bot, l_right, radius)[0])
arcs.append(part.fillet_lines(l_right, l_top, radius)[0])
model.set_ediv(arcs, 10)

part.chunk()
model.plot_geometry(pnum=False, lnum=False) # view the geometry

model.set_etype('plstress', part, 0.01)
model.set_eshape('quad', 2)
model.mesh(0.7, 'gmsh')
model.plot_elements()
