#!/usr/bin/env python3
import pycalculix as pyc

model_name = 'pinned-plate'
model = pyc.FeaModel(model_name)
model.set_units('in')

# pin locations
pin1 = [0, 0]
pin2 = [pin1[0], 4]
pin3 = [4, 8]

# dimentions for drawing parts
pinhole_rad = 0.25
pin_rad = pinhole_rad - .1
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
hole_arcs = []
for pin_point in [pin1, pin2, pin3]:
    hole = part.draw_hole(pin_point[0], pin_point[1], pinhole_rad)
    hole_arcs += hole

# make pins
pins = []
pin_arcs = []
for pin_point in [pin1, pin2, pin3]:
    pin = pyc.Part(model)
    arcs = pin.draw_circle(pin_point[0], pin_point[1], pin_rad)
    pin.chunk(exclude_convex=False)
    pins.append(pin)
    pin_arcs += arcs

all_arcs = fillet_arcs + hole_arcs + pin_arcs
model.set_ediv(all_arcs, 10)

disp = False
#pyc.feamodel.GEOM_CMAP = 'Set2'
pyc.feamodel.GEOM_CMAP = 'Paired'
model.plot_parts(model_name+'_prechunk_parts')
model.plot_areas(model_name+'_prechunk_areas', label=False)
model.plot_lines(model_name+'_prechunk_lines', display=disp)
model.plot_points(model_name+'_prechunk_areas', display=disp)
model.plot_geometry(model_name+'_prechunk_geom', display=disp)
part.chunk('ext')
model.plot_areas(model_name+'_areas', label=False)
model.plot_lines(label=False, display=disp)

"""
model.set_eshape('tri', 2)
model.set_etype('plstress', part, 0.1)
model.set_etype('plstress', pins, 0.1)
model.mesh(0.5, 'gmsh')
model.plot_elements()
"""


"""
import pickle
import sys
sys.setrecursionlimit(100000)
pickle.dump(model, open('test.pickle', 'wb'))
print('pickled')
del model
model = pickle.load(open('test.pickle', 'rb'))
print('loaded pickle')
model.plot_geometry()
"""