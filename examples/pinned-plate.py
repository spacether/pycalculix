#!/usr/bin/env python3
import pycalculix as pyc

# Model of a pinned plate with 3 pins
# one will be force, the others will be fixed
proj_name = 'pinned-plate'
model = pyc.FeaModel(proj_name)
model.set_units('in')

# pin locations
pin1 = [0, 0]
pin2 = [pin1[0], 4]
pin3 = [4, 8]

# dimensions for drawing parts
pinhole_rad = 0.25
pin_rad = pinhole_rad - .015
width_plate = 1
near_hole_th = 0.5*width_plate
left_y = pin1[1] - near_hole_th
left_x = pin1[0] - near_hole_th
pinhole_dist = pin2[1] - pin1[1]

# make main part
part = pyc.Part(model)
part.goto(left_x, left_y)
part.draw_line_ax(pinhole_dist + width_plate)
part.draw_line_rad(width_plate)
part.draw_line_ax(-width_plate)
part.draw_line_ax(-(pinhole_dist-width_plate))
part.draw_line_rad(pin3[0] - pin1[0] - width_plate)
part.draw_line_ax(pin3[1] - pin1[0])
part.draw_line_rad(width_plate)
part.draw_line_ax(-width_plate)
part.draw_line_ax(-(pin3[1] - pin1[1]))
part.draw_line_to(left_x, left_y)
fillet_arcs = part.fillet_all(near_hole_th)

# make pinholes
holes = []
hole_arcs = []
for pin_point in [pin1, pin2, pin3]:
    hole = part.draw_hole(pin_point[0], pin_point[1], pinhole_rad)
    hole_arcs += hole
    holes.append(hole)

# make pins
pin_parts = []
pins = []
pin_arcs = []
for pin_point in [pin1, pin2, pin3]:
    pin = pyc.Part(model)
    arcs = pin.draw_circle(pin_point[0], pin_point[1], pin_rad)
    #pin.chunk(exclude_convex=False)
    pins.append(arcs)
    pin_parts.append(pin)
    pin_arcs += arcs

# set divisions on pin holes and pins
all_arcs = fillet_arcs + hole_arcs + pin_arcs
model.set_ediv(all_arcs, 10)

# plot model
model.plot_areas(proj_name+'_prechunk_areas', label=False)
part.chunk('ext')
model.plot_areas(proj_name+'_areas', label=False)
model.plot_parts(proj_name+'_parts')
model.plot_points(proj_name+'_points')
model.plot_lines(proj_name+'_points', label=False)

# set loads and constraints
pin_vert_pts = model.get_items(['P40','P42', 'P45', 'P47'])
pin_horiz_pts = model.get_items(['P38','P41', 'P43', 'P46'])
top_pin_vert_pts = model.get_items(['P50', 'P52'])
top_pin_pt = model.get_item('P50')

model.set_constr('fix',pin_horiz_pts,'x')
model.set_constr('fix',pin_vert_pts,'y')
model.set_constr('fix',top_pin_vert_pts,'y')
model.set_load('force',top_pin_pt,-500,'x')
model.set_gravity(9.81, [part] + pin_parts)

# set part material
mat = pyc.Material('steel')
youngs = 210*(10**9)
mat.set_mech_props(7800, youngs, 0.3)
model.set_matl(mat, pin_parts)
model.set_matl(mat, part)

# set contact
factor = 5 # can be between 5 and 50
kval = youngs*factor
for (pin, hole) in zip(pins, holes):
    model.set_contact_linear(pin, hole, kval, True)

# mesh the model
model.set_eshape('tri', 2)
model.set_etype('plstress', part, 0.1)
model.set_etype('plstress', pin_parts, 0.1)
model.mesh(0.5, 'gmsh')
model.plot_elements(proj_name+'_elem')   # plot part elements
model.plot_constraints(proj_name+'_constr')

# make and solve the model
prob = pyc.Problem(model, 'struct')
prob.solve()

# Plot results
fields = 'Sx,Sy,S1,S2,S3,Seqv,ux,uy,utot'    # store the fields to plot
fields = fields.split(',')
for field in fields:
    fname = proj_name+'_'+field
    prob.rfile.nplot(field, fname, display=False)

model.view.select(part)
model.view.allsel_under('parts')
for field in fields:
    fname = proj_name+'_PART_'+field
    prob.rfile.nplot(field, fname, display=False)
