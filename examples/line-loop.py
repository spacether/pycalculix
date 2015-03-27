#!/usr/bin/env python3
import pycalculix as pyc

# CLOCKWISE square
coords = [[0,0],
          [3,0],
          [3,3],
          [0,3]]
coords.append(coords[0])

lines = []
for ind in range(1,len(coords)):
    start = coords[ind-1]
    end = coords[ind]
    start = pyc.geometry.Point(start[0], start[1])
    end = pyc.geometry.Point(end[0], end[1])
    lines.append(pyc.geometry.Line(start, end))
print('-----------------------------------------')
print('Clockwise square')
loop = pyc.geometry.LineLoop(lines)
print('area: %f' % loop.area)
print('center: %s' % loop.center)
print('ccw: %s' % loop.ccw)

# CLOCKWISE square
coords = [[1,1],
          [2,1],
          [2,2],
          [1,2]]
coords.append(coords[0])

lines = []
for ind in range(1,len(coords)):
    start = coords[ind-1]
    end = coords[ind]
    start = pyc.geometry.Point(start[0], start[1])
    end = pyc.geometry.Point(end[0], end[1])
    lines.append(pyc.geometry.Line(start, end))
print('-----------------------------------------')
print('Clockwise square')
small_loop = pyc.geometry.LineLoop(lines)
print('area: %f' % small_loop.area)
print('center: %s' % small_loop.center)
print('ccw: %s' % small_loop.ccw)

pts = set()
for sline in small_loop:
    pts.update(sline.points)
print ('small_loop points:')
for point in pts:
    print(point)
print('Checking if small is inside big')
patch = loop.get_patch()
for point in pts:
    contains = patch.contains_point((point.x, point.y))
    print('Checking point %s' % point)
    print('Contains: %s' % contains) 
inside = small_loop.inside(loop)
print('Small is inside big: %s' % inside)

"""
# CLOCKWISE
# radial, ax
coords = [[0,0],
          [2,0],
          [2,1],
          [[1,2],[1,1]],
          [0,2]]
coords.append(coords[0])

lines = []
for ind in range(1,len(coords)):
    start = coords[ind-1]
    end = coords[ind]
    actr = None
    if isinstance(start[0], list):
        # draw line from arc
        start = start[0]
    if isinstance(end[0], list):
        # draw arc
        actr = end[1]
        actr = pyc.geometry.Point(actr[0], actr[1])
        end = end[0]
    start = pyc.geometry.Point(start[0], start[1])
    end = pyc.geometry.Point(end[0], end[1])
    if actr == None:
        # line
        lines.append(pyc.geometry.Line(start, end))
    else:
        lines.append(pyc.geometry.Arc(start, end, actr))
print('-----------------------------------------')
print('Clockwise square with filleted corner')
loop = pyc.geometry.LineLoop(lines)
print('area: %f' % loop.area)
print('center: %s' % loop.center)

# COUNTER-CLOCKWISE
# radial, ax
coords = [[0,0],
          [0,2],
          [1,2],
          [[2,1],[1,1]],
          [2,0]]
coords.append(coords[0])

lines = []
for ind in range(1,len(coords)):
    start = coords[ind-1]
    end = coords[ind]
    actr = None
    if isinstance(start[0], list):
        # draw line from arc
        start = start[0]
    if isinstance(end[0], list):
        # draw arc
        actr = end[1]
        actr = pyc.geometry.Point(actr[0], actr[1])
        end = end[0]
    start = pyc.geometry.Point(start[0], start[1])
    end = pyc.geometry.Point(end[0], end[1])
    if actr == None:
        # line
        lines.append(pyc.geometry.Line(start, end))
    else:
        lines.append(pyc.geometry.Arc(start, end, actr))
print('-----------------------------------------')
print('Counter-Clockwise square with filleted corner')
loop = pyc.geometry.LineLoop(lines)
print('area: %f' % loop.area)
print('center: %s' % loop.center)


# CLOCKWISE
# radial, ax
coords = [[0,1],
          [[1,0],[0,0]],
          [[0,-1],[0,0]],
          [[-1,0],[0,0]],
          [[0,1],[0,0]]]

lines = []
for ind in range(1,len(coords)):
    start = coords[ind-1]
    end = coords[ind]
    actr = None
    if isinstance(start[0], list):
        # draw line from arc
        start = start[0]
    if isinstance(end[0], list):
        # draw arc
        actr = end[1]
        actr = pyc.geometry.Point(actr[0], actr[1])
        end = end[0]
    start = pyc.geometry.Point(start[0], start[1])
    end = pyc.geometry.Point(end[0], end[1])
    if actr == None:
        # line
        lines.append(pyc.geometry.Line(start, end))
    else:
        lines.append(pyc.geometry.Arc(start, end, actr))
print('-----------------------------------------')
print('Counter-Clockwise circle')
loop = pyc.geometry.LineLoop(lines)
print('area: %f' % loop.area)
print('center: %s' % loop.center)
print('ccw: %r' % loop.ccw)

# CLOCKWISE
# radial, ax
coords = [[0,0],
          [2,0],
          [2,1],
          [[1,2],[2,2]],
          [0,2]]
coords.append(coords[0])

lines = []
for ind in range(1,len(coords)):
    start = coords[ind-1]
    end = coords[ind]
    actr = None
    if isinstance(start[0], list):
        # draw line from arc
        start = start[0]
    if isinstance(end[0], list):
        # draw arc
        actr = end[1]
        actr = pyc.geometry.Point(actr[0], actr[1])
        end = end[0]
    start = pyc.geometry.Point(start[0], start[1])
    end = pyc.geometry.Point(end[0], end[1])
    if actr == None:
        # line
        lines.append(pyc.geometry.Line(start, end))
    else:
        lines.append(pyc.geometry.Arc(start, end, actr))
print('-----------------------------------------')
print('Clockwise square with concave corner')
loop = pyc.geometry.LineLoop(lines)
print('area: %f' % loop.area)
print('ccw: %r' % loop.ccw)
print('center: %s' % loop.center)
loop.reverse()
print('area: %f' % loop.area)
print('ccw: %r' % loop.ccw)
print('center: %s' % loop.center)
"""