from pycalculix import FeaModel
import matplotlib.pyplot as plt
import math

# We'll be modeling a rotating jet engine part
# Problem constants
proj_name = 'example3_compr_rotor'

# problem + geometry constants
rpm = 1000     # rotor speed in rpm
rad_inner = 1
rad_flowpath = 3
disk_width = 1
disk_ht = 0.5
disk_angle = 45 # from horizontal, internal to disk
web_width = 0.2
arm_length = .5
arm_th = 0.05
airfoil_width = 0.4
airfoil_ht = 1
rad_disk_top = 0.1
rad_disk_bot = 0.1
rad_web_top = 0.1
rad_web_bot = 0.1

# derived dimensions
diskramp_ax = (disk_width - web_width)/2
diskramp_rad = math.tan(math.radians(disk_angle))*diskramp_ax
web_ht = rad_flowpath - arm_th - diskramp_rad - rad_inner
flowpath_ax = (web_width + 2*arm_length - airfoil_width) / 2.0

# these are deltas
disk_loop = [[disk_ht,0],[diskramp_rad,diskramp_ax],[web_ht,0,],[0,-arm_length],
             [arm_th,0],[0,flowpath_ax],[0,airfoil_width],[0,flowpath_ax],
             [-arm_th,0],[0,-arm_length],[-web_ht,0],[-diskramp_rad,diskramp_ax],
             [-disk_ht,0],[0,-disk_width]]

airfoil_loop = [[airfoil_ht,0],[0,airfoil_width],[-airfoil_ht,0],[0,-airfoil_width]]

# make model
a = FeaModel(proj_name)
a.set_units('m')    # this sets dist units to meters, labels our consistent units

# make part
b = a.PartMaker()
b.goto(rad_inner,0)
lines = []
for [drad, dax] in disk_loop:
    [L, p1, p2] = b.draw_line_delta(drad, dax)
    lines.append(L)
# go to where we'll draw the airfoil
mypt = lines[6].pt(0) # first point of line 8
b.goto(mypt.x,mypt.y)
for [drad, dax] in airfoil_loop:
    b.draw_line_delta(drad, dax)

# fillets [[line1, line2, fillet_rad],[
fillet_list = [[13,0,rad_disk_bot],[12,13,rad_disk_bot],
               [0,1,rad_disk_top],[11,12,rad_disk_top],
               [1,2,rad_web_bot],[10,11,rad_web_bot],
               [2,3,rad_web_top],[9,10,rad_web_top]]

for [i1, i2, rad] in fillet_list:
    b.fillet_lines( lines[i1], lines[i2], rad )

a.plot_geometry(proj_name+'_geom') # view the geometry, points, lines, and areas

# set loads and constraints
a.set_rpm(10000, b)
a.set_constr('fix',b.left,'y')

# set part material
mat = a.MatlMaker('nickel_alloy')
mat.set_mech_props(8220, 208*10**9, 0.3)
a.set_matl(mat, b)

# mesh
a.set_eshape('quad', 2)
a.set_etype(b, 'axisym')
a.set_etype('A1', 'plstress', 0.1)
a.get_item('L15').set_ediv(8)
a.get_item('L6').set_ediv(8)
a.get_item('L2').set_ediv(24)
a.get_item('L10').set_ediv(24)
a.mesh(1.0, 'gmsh') # mesh with 1.0 fineness, smaller is finer
a.plot_elements(proj_name+'_elem')   # plot the part elements
a.plot_constraints(proj_name+'_constr')

# make and solve the model
mod = a.ModelMaker(b, 'struct')
mod.solve()

# view the results
disp = False    # turn off display plotting
g = 0.0         # turn off adding results displacement to the model
fields = 'Sx,Sy,Sz,S1,S2,S3,Seqv,ux,uy,uz,utot'    # store the fields to plot
fields = fields.split(',')
for field in fields:
    fname = proj_name+'_'+field
    mod.rfile.nplot(field, fname, display=disp, gmult=g)
