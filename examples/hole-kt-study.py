#!/usr/bin/env python3
import math
import numpy as np
import sys

import pycalculix as pyc
import matplotlib.pyplot as plt

# set whether or not to show gui plots
show_gui = True
if '-nogui' in sys.argv:
    show_gui = False
# set element shape
eshape = 'quad'
if '-tri' in sys.argv:
    eshape = 'tri'

# Stress and geometry constants
stress_val = 1000
diam = 1.0
thickness = 0.01

# Make a list of geometry ratios, diam_hole/width_plate
ratios = np.arange(0,.5,.05)
ratios[0] = .001

# store results
(ktg_res, ktg_pet, err) = ([],[],[])

def kt_peterson(ratio):
    # returns peterson kt for a given ratio, kt is Ktg
    res = .284 + (2.0/(1-ratio)) - 0.600*(1-ratio) + 1.32*(1-ratio)**2
    return res

# loop through ratios doing a stress run each time, storing the results
for ratio in ratios:
    width = diam/ratio

    # part geometry dimensions
    print('D=%f, H=%f, D/H=%f' % (diam, width, diam/width))
    top = width/2   # model width
    right = top*2  # model length
    rad = diam/2.0   #hole radius
    bot = top - rad
    left = right - rad

    # vertical hole in plate model, make model
    model_name = 'hole-kt-study'
    model = pyc.FeaModel(model_name)
    model.set_units('m')    # this sets dist units to meters

    # make part, coordinates are x, y = radial, axial
    part = pyc.Part(model)
    part.goto(0.0,rad)
    part.draw_arc(rad, 0.0, 0.0, 0.0)
    part.draw_line_rad(left)
    part.draw_line_ax(top)
    part.draw_line_rad(-right*.5)
    part.draw_line_rad(-right*.5) #this point lets us chunks our area
    part.draw_line_ax(-bot)
    # part.plot_geometry('hole_kt_prechunk', display=disp)
    part.chunk()
    model.plot_geometry(model_name+'_chunked', display=False)

    # set loads and constraints
    model.set_load('press',part.top,-1*stress_val)
    model.set_constr('fix',part.left,'y')
    model.set_constr('fix',part.bottom,'x')

    # set part material
    mat = pyc.Material('steel')
    mat.set_mech_props(7800, 210000, 0.3)
    model.set_matl(mat, part)

    # set the element type, line division, and mesh the database
    ediv = 19
    model.set_ediv('L0',ediv) # sets # of elements on the arc

    model.set_eshape(eshape, 2)
    model.set_etype('plstress', part, thickness)
    model.mesh(1.0, 'gmsh') # mesh with 1.0 fineness, smaller is finer
    model.plot_elements('%s_elem_%.3f' % (model_name, ratio),
                        display=False)
    model.plot_pressures('%s_press' % (model_name), display=False)

    # make model and solve it
    prob = pyc.Problem(model, 'struct')
    prob.solve()

    # query results and store them
    sx = prob.rfile.get_nmax('Sx')
    kt_fea = sx/stress_val
    ktg_res.append(kt_fea)
    ktg_pet.append(kt_peterson(ratio))
    error = 100*(kt_fea/kt_peterson(ratio) - 1)
    err.append(error)
    print('For ratio %3f, Kt_g = %3.2f' % (ratio, kt_fea))

# plot results
fig, ax = plt.subplots()
plt.plot(ratios, ktg_res, color='b', label='Ktg_FEA', marker='.')
plt.plot(ratios, ktg_pet, color='r', label='Ktg_Peterson', marker='.')
plt.grid()
plt.legend(loc='lower right')
plt.title('Tension Hole in Plate Stress Concentration Factor, Ktg')
plt.xlabel('D/h')
plt.ylabel('Ktg')
pyc.base_classes.plot_finish(plt, fname=model_name+'_kts',
                             display=show_gui)

# plot error
fig, ax = plt.subplots()
plt.plot(ratios, err, color='g', label='Error', marker='.')
plt.grid()
plt.legend(loc='lower right')
plt.title('Tension Hole in Plate Ktg Error, FEA vs Peterson')
plt.xlabel('D/h')
plt.ylabel('Error (%)')
pyc.base_classes.plot_finish(plt, fname=model_name+'_error',
                             display=show_gui)
