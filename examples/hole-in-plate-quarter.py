#!/usr/bin/env python3
import sys

import pycalculix as pyc

# Vertical hole in plate model, make model
proj_name = 'hole-in-plate-quarter'
model = pyc.FeaModel(proj_name)
model.set_units('m') # this sets dist units to meters

# the below boolean sets whether or not to show gui plots
# testing passes in -nogui
show_gui = True
if len(sys.argv) == 2 and sys.argv[-1] == '-nogui':
    show_gui = False

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
model.plot_geometry(proj_name+'_A0', display=show_gui)
part.chunk() # cut the part into area pieces so CGX can mesh it
model.plot_geometry(proj_name+'_A0chunked', display=show_gui)

# set loads and constraints
model.set_load('press',part.top,-1000)
model.set_constr('fix',part.left,'y')
model.set_constr('fix',part.bottom,'x')

# set part material
mat = pyc.Material('steel')
mat.set_mech_props(7800, 210*(10**9), 0.3)
model.set_matl(mat, part)

# set the element type and mesh database
model.set_eshape('quad', 2)
model.set_etype('plstress', part, 0.1)
model.set_ediv('L0', 20) # set element divisions
model.mesh(1.0, 'gmsh') # mesh 1.0 fineness, smaller is finer
model.plot_elements(proj_name+'_elem', display=show_gui)
model.plot_pressures(proj_name+'_press', display=show_gui)
model.plot_constraints(proj_name+'_constr', display=show_gui)

# make and solve the model
prob = pyc.Problem(model, 'struct')
prob.solve()

# view and query results
sx = prob.rfile.get_nmax('Sx')
print('Sx_max: %f' % sx)
[fx, fy, fz] = prob.rfile.get_fsum(model.get_item('L5'))
print('Reaction forces (fx,fy,fz) = (%12.10f, %12.10f, %12.10f)' %
      (fx, fy, fz))

# Plot results
# store the fields to plot
fields = 'Sx,Sy,S1,S2,S3,Seqv,ux,uy,utot,ex'
fields = fields.split(',')
for field in fields:
    fname = proj_name+'_'+field
    prob.rfile.nplot(field, fname, display=False)
