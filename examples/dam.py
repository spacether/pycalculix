#!/usr/bin/env python3
import math
import sys

import pycalculix as pyc

# We'll be modeling a masonry gravity dam, the Beetaloo dam in Australia
# Problem constants
proj_name = 'dam'
grav = 9.81 # m/s^2
dens_water = 1000 #kg/m^3
press_atm = 101325 # Pascals = N/m^2 = kg/(m-s^2)
dam_ht_ft = 115
water_ht_ft = 109
length_dam_ft = 580

# the below boolean sets whether or not to show gui plots
# testing passes in -nogui
show_gui = True
if len(sys.argv) == 2 and sys.argv[-1] == '-nogui':
    show_gui = False

# derived dims
water_ax = (water_ht_ft-22)/math.tan(math.radians(89)) + 12
top_ax = (dam_ht_ft-22)/math.tan(math.radians(89)) + 12
thickness = length_dam_ft*0.3048 # conversion to metric

# derived dimensions, tan = o/a --> o/tan = a
pts_ft_water = [[8,0],[22,12],[water_ht_ft,water_ax]]
pts_ft_air = [[dam_ht_ft,top_ax], [dam_ht_ft, top_ax + 14]]
water_ht_m = water_ht_ft*0.3048

# make model
model = pyc.FeaModel(proj_name)
# this sets dist units to meters, labels our consistent units
model.set_units('m')

# make part, coordinates are x, y = radial, axial
part = pyc.Part(model)
part.goto(0.0,0.0)
water_lines = []
air_lines = []
for [x,y] in pts_ft_water:
    [x,y] = [x*0.3048,y*0.3048] # conversion to metric
    [L1,p1,p2] = part.draw_line_to(x, y)
    water_lines.append(L1)
for [x,y] in pts_ft_air:
    [x,y] = [x*0.3048,y*0.3048] # conversion to metric
    [L1,p1,p2] = part.draw_line_to(x, y)
    air_lines.append(L1)
# make the two arcs
pts_ft_arcs = [ [[22,73],[146,208]], [[14,98],[41,93]] ]
for [[x,y],[xc,yc]] in pts_ft_arcs:
    [x,y] = [x*0.3048,y*0.3048]
    [xc,yc] = [xc*0.3048,yc*0.3048]
    [L1,p1,p2] = part.draw_arc(x,y,xc,yc)
    air_lines.append(L1)
# make the last 3 lines after the arcs
pts_ft_other = [ [2,110], [0, 110] ]
for [x,y] in pts_ft_other:
    [x,y] = [x*0.3048,y*0.3048] # conversion to metric
    [L1,p1,p2] = part.draw_line_to(x, y)
    air_lines.append(L1)
part.draw_line_to(0, 0)
# view the points, lines, and areas
model.plot_geometry(proj_name+'_geom', display=show_gui)

# set part material
mat = pyc.Material('concrete')
mat.set_mech_props(2300, 30000*(10**6), 0.2)
model.set_matl(mat, part)

# set the element type, line division, and mesh the database
model.set_eshape('quad', 2)
model.set_etype('plstrain', part, thickness)
model.get_item('L8').set_ediv(2)
# mesh with 1.0 or less fineness, smaller is finer
model.mesh(0.5, 'gmsh')
# plot the part elements
model.plot_elements(proj_name+'_elem', display=show_gui)

# set loads and constraints
model.set_load('press', air_lines, press_atm)
model.set_fluid_press(water_lines, dens_water, grav, water_ht_m,
                      press_atm)
model.set_gravity(grav, part)
model.set_constr('fix', part.bottom, 'x')
model.set_constr('fix', part.bottom, 'y')
model.plot_pressures(proj_name+'_press_1', display=show_gui)
model.plot_constraints(proj_name+'_constr', display=show_gui)

# make model and solve it
prob = pyc.Problem(model, 'struct')
prob.solve()

# query results and store them
disp = False  # turn off display plotting
# store the fields to write
fields = 'Seqv,Sx,Sy,Sz,S1,S2,S3,ux,uy,utot'
fields = fields.split(',')
for field in fields:
    fname = proj_name+'_'+field
    prob.rfile.nplot(field, fname, display=disp)
smax = prob.rfile.get_nmax('Seqv')
[fx, fy, fz] = prob.rfile.get_fsum(model.get_item('L9'))
print('Seqv_max= %3.2f' % (smax))
print('Reaction forces (fx,fy,fz) = (%12.10f, %12.10f, %12.10f)' % (fx, fy, fz))
