#!/usr/bin/env python3
import pycalculix as pyc

# Model of a pipe being crushed, quarter-symmetrym plane strain
proj_name = 'pipe-crush-elastic'
model = pyc.FeaModel(proj_name)
model.set_units('m') # this sets dist units to meters

# Define variables we'll use to draw part geometry
rad_outer = 0.1143 # pipe outer radius
pipe_wall_th = .00887 # wall thickness
rad_inner = rad_outer - pipe_wall_th
pipe_length = 10*rad_outer
e_size = pipe_wall_th/4.0

# Draw pipe, x, y = radial, axial
pipe = pyc.Part(model)
pipe.goto(0.0, rad_inner)
arc_inner = pipe.draw_arc_angle(90, 0.0, 0.0)[0]
pipe.draw_line_to(rad_outer, 0.0)
arc_outer = pipe.draw_arc_angle(-90, 0.0, 0.0)[0]
pipe.draw_line_to(0.0, rad_inner)
num_outer = arc_outer.length()/e_size
num_inner = arc_inner.length()/e_size
num_arc_eles = int(round(0.5*(num_outer + num_inner),0))

# Draw plate, x, y = radial, axial
plate = pyc.Part(model)
plate.goto(rad_outer, 0.0)
plate.draw_line_rad(pipe_wall_th)
plate.draw_line_ax(rad_outer)
plate.draw_line_rad(-pipe_wall_th)
plate_bottom = plate.draw_line_to(rad_outer, 0.0)[0]
num_plate_eles = int(round(plate_bottom.length()/e_size,0))

# view model
model.plot_geometry(proj_name+'_geometry')
model.plot_parts(proj_name+'_parts')
model.plot_areas(proj_name+'_areas')
model.plot_lines(proj_name+'_lines')

# set loads and constraints
model.set_constr('fix',pipe.left,'y')
model.set_constr('fix',pipe.bottom,'x')
model.set_constr('fix',plate.left,'y')
ux_total = -.050
percent_displ = 50
model.set_constr('displ',plate.top,'x',ux_total*percent_displ/100.0)

# set part material
mat = pyc.Material('steel')
youngs = 210*(10**9)
mat.set_mech_props(7800, youngs, 0.3)
model.set_matl(mat, pipe)
model.set_matl(mat, plate)

# set contact
factor = 5 # can be between 5 and 50
kval = youngs*factor
model.set_contact_linear(plate_bottom, arc_outer, kval)

# set the element type and mesh database
model.set_eshape('quad', 2)
model.set_etype('plstrain', pipe, pipe_length)
model.set_etype('plstrain', plate, pipe_length)
model.set_ediv(['L1','L3', 'L4', 'L6'], 4) # set element divisions
model.set_ediv(['L0','L2'], num_arc_eles) # set element divisions
model.set_ediv(['L5','L7'], num_plate_eles)
model.mesh(1.0, 'gmsh') # mesh 1.0 fineness, smaller is finer
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
