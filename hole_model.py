from pycalculix import FeaModel

# Vertical hole in plate model, make model
proj_name = 'hole_model'
a = FeaModel(proj_name)
a.set_units('m')    # this sets dist units to meters, labels our consistent units

# Define variables we'll use to draw part geometry
diam = 2.0 # hole diam
ratio = 0.45
width = diam/ratio   #plate width
print('D=%f, H=%f, D/H=%f' % (diam, width, diam/width))
length = 2*width  #plate length
rad = diam/2   #hole radius
vdist = (length - 2*rad)/2  #derived dimension
adist = width/2             #derived dimension

# Draw part geometry, you must draw the part CLOCKWISE
# coordinates are x, y = radial, axial
b = a.PartMaker()
b.goto(0.0,rad)
b.draw_arc(rad, 0.0, 0.0, 0.0)
b.draw_line_rad(vdist)
b.draw_line_ax(adist)
b.draw_line_rad(-length/4.0)
b.draw_line_rad(-length/4.0)
b.draw_line_ax(-(adist-rad))
a.plot_geometry(proj_name+'_prechunk') # view the geometry

# Cut the part into easier to mesh areas
b.chunk() # cut the part into area pieces so CGX can mesh it
a.plot_geometry(proj_name+'_chunked') # view the geometry

# set loads and constraints
a.set_load('press',b.top,-1000)
a.set_constr('fix',b.left,'y')
a.set_constr('fix',b.bottom,'x')

# set part material
mat = a.MatlMaker('steel')
mat.set_mech_props(7800, 210000, 0.3)
a.set_matl(mat, b)

# set the element type and mesh database
a.set_eshape('quad', 2)
a.set_etype(b, 'plstress', 0.1)
a.get_item('L0').set_ediv(20) # set element divisions
a.mesh(1.0, 'gmsh') # mesh 1.0 fineness, smaller is finer
a.plot_elements(proj_name+'_elem')   # plot part elements
a.plot_pressures(proj_name+'_press')
a.plot_constraints(proj_name+'_constr')

# make and solve the model
mod = a.ModelMaker(b, 'struct')
mod.solve()

# view and query results
sx = mod.rfile.get_nmax('Sx')
print('Sx_max: %f' % (sx))
[fx, fy, fz] = mod.rfile.get_fsum(b.get_item('L5'))
print('Reaction forces (fx,fy,fz) = (%12.10f, %12.10f, %12.10f)' % (fx, fy, fz)) 

# Plot results
disp = False
mod.rfile.nplot('Sx', proj_name+'_Sx', display=disp)
mod.rfile.nplot('Sy', proj_name+'_Sy', display=disp)
mod.rfile.nplot('S1', proj_name+'_S1', display=disp)
mod.rfile.nplot('S2', proj_name+'_S2', display=disp)
mod.rfile.nplot('S3', proj_name+'_S3', display=disp)
mod.rfile.nplot('Seqv', proj_name+'_Seqv', display=disp)
mod.rfile.nplot('ux', proj_name+'_ux', display=disp)
mod.rfile.nplot('uy', proj_name+'_uy', display=disp)
mod.rfile.nplot('utot', proj_name+'_utot', display=disp)

