import pycalculix as pyc

# Vertical hole in plate model, make model
proj_name = 'hole-in-plate-full'
model = pyc.FeaModel(proj_name)
model.set_units('m') # this sets dist units to meters

# Define variables we'll use to draw part geometry
diam = 2.0 # hole diam
ratio = 0.45
width = diam/ratio   #plate width
print('D=%f, H=%f, D/H=%f' % (diam, width, diam/width))
length = 2*width  #plate length

# Draw part geometry, you must draw the part CLOCKWISE, x, y = radial, axial
part = pyc.Part(model)
part.goto(length*0.5, -width*0.5)
part.draw_line_ax(width)
part.draw_line_rad(-length)
part.draw_line_ax(-width)
part.draw_line_rad(length)
hole_lines = part.draw_hole(0, 0, diam*0.5, filled=False)
model.set_ediv(hole_lines, 10)
part.chunk()
model.plot_geometry(proj_name+'_geom') # view the geometry

# set loads and constraints
pressure = -1000
model.set_load('press', part.top, pressure)
model.set_load('press', part.bottom, pressure)
model.set_constr('fix', ['P6', 'P8'], 'y')
model.set_constr('fix', ['P4', 'P7'], 'x')

# set part material
mat = pyc.Material('steel')
mat.set_mech_props(7800, 210*(10**9), 0.3)
model.set_matl(mat, part)

# set the element type and mesh database
model.set_eshape('quad', 2)
model.set_etype('plstress', part, 0.1)
model.mesh(1.0, 'gmsh') # mesh 1.0 fineness, smaller is finer
model.plot_elements(proj_name+'_elem')   # plot part elements
model.plot_pressures(proj_name+'_press')
model.plot_constraints(proj_name+'_constr')

# make and solve the model
prob = pyc.Problem(model, 'struct')
prob.solve()

# view and query results
sx = prob.rfile.get_nmax('Sx')
print('Sx_max: %f' % sx)

# Plot results
fields = 'Sx,Sy,S1,S2,S3,Seqv,ux,uy,utot,ex'    # store the fields to plot
fields = fields.split(',')
for field in fields:
    fname = proj_name+'_'+field
    prob.rfile.nplot(field, fname, display=False)
