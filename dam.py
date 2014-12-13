from pycalculix import FeaModel, frange
import matplotlib.pyplot as plt
import math

# We'll be modeling a masonry gravity dam, the Beetaloo dam in Australia

# Problem constants
proj_name = 'dam'
grav = 9.81 # m/s^2
dens_water = 1000 #kg/m^3
press_atm = 101325 # Pascals = N/m^2 = kg/(m-s^2)
dam_ht_ft = 115
water_ht_ft = 109
length_dam_ft = 580

# derived dims
water_ax = (water_ht_ft-22)/math.tan(math.radians(89)) + 12
top_ax = (dam_ht_ft-22)/math.tan(math.radians(89)) + 12
thickness = length_dam_ft*0.3048 # conversion to metric

# derived dimensions, tan = o/a --> o/tan = a
pts_ft_water = [[8,0],[22,12],[water_ht_ft,water_ax]]
pts_ft_air = [[dam_ht_ft,top_ax], [dam_ht_ft, top_ax + 14]]
water_ht_m = water_ht_ft*0.3048

# make model
a = FeaModel(proj_name)
a.set_units('m')    # this sets dist units to meters, labels our consistent units

# make part, coordinates are x, y = radial, axial
b = a.PartMaker()
b.goto(0.0,0.0)
water_lines = []
air_lines = []
for [x,y] in pts_ft_water:
    [x,y] = [x*0.3048,y*0.3048] # conversion to metric
    [L1,p1,p2] = b.draw_line_to(x, y)
    water_lines.append(L1)
for [x,y] in pts_ft_air:
    [x,y] = [x*0.3048,y*0.3048] # conversion to metric
    [L1,p1,p2] = b.draw_line_to(x, y)
    air_lines.append(L1)   
# make the two arcs
pts_ft_arcs = [ [[22,73],[146,208]], [[14,98],[41,93]] ]
for [[x,y],[xc,yc]] in pts_ft_arcs:
    [x,y] = [x*0.3048,y*0.3048]
    [xc,yc] = [xc*0.3048,yc*0.3048]
    [L1,p1,p2] = b.draw_arc(x,y,xc,yc)
    air_lines.append(L1)
# make the last 3 lines after the arcs
pts_ft_other = [ [2,110], [0, 110] ]
for [x,y] in pts_ft_other:
    [x,y] = [x*0.3048,y*0.3048] # conversion to metric
    [L1,p1,p2] = b.draw_line_to(x, y)
    air_lines.append(L1)
b.draw_line_to(0, 0)
a.plot_geometry(proj_name+'_geom', b) # view the geometry, points, lines, and areas

# set part material
mat = a.MatlMaker('concrete')
mat.set_mech_props(2300, 30000*(10**6), 0.2)
a.set_matl(mat, b)

# set the element type, line division, and mesh the database
a.set_eshape('quad', 2)
a.set_etype(b, 'plstrain', thickness)
b.get_item('L8').set_ediv(2)
a.mesh(0.5, 'gmsh')               # mesh with 1.0 fineness, smaller is finer
a.plot_elements(proj_name+'_elem', b)   # plot the part elements

# set loads and constraints
a.set_load('press', air_lines, press_atm)
a.set_fluid_press(water_lines, dens_water, grav, water_ht_m, press_atm)
a.set_gravity(grav, b)
a.set_constr('fix', b.bottom, 'x')
a.set_constr('fix', b.bottom, 'y')
a.plot_pressures(proj_name+'_press', b)

'''
a.set_time(2.0)
a.set_load('press', water_lines, press_atm)
'''

# make model and solve it
mod = a.ModelMaker(b, 'struct')
mod.solve()

# query results and store them
disp = False    # turn off display plotting
mod.rfile.nplot('Seqv', proj_name+'_Seqv', display=disp)
mod.rfile.nplot('Sx', proj_name+'_Sx', display=disp)
mod.rfile.nplot('Sy', proj_name+'_Sy', display=disp)
mod.rfile.nplot('Sz', proj_name+'_Sz', display=disp)
mod.rfile.nplot('S1', proj_name+'_S1', display=disp)
mod.rfile.nplot('S2', proj_name+'_S2', display=disp)
mod.rfile.nplot('S3', proj_name+'_S3', display=disp)
mod.rfile.nplot('ux', proj_name+'_ux', display=disp)
mod.rfile.nplot('uy', proj_name+'_uy', display=disp)
mod.rfile.nplot('utot', proj_name+'_utot', display=disp)

mod.rfile.nplot('ux', proj_name+'_ux_levels', display=disp)
mod.rfile.nplot('ux', proj_name+'_ux_gradient', display=disp, gradient=True)

'''
smax = mod.rfile.get_nmax('Seqv')
[fx, fy, fz] = mod.rfile.get_fsum(a.get_item('L9'))
print('Seqv_max= %3.2f' % (smax))
print('Reaction forces (fx,fy,fz) = (%12.10f, %12.10f, %12.10f)' % (fx, fy, fz)) 
'''