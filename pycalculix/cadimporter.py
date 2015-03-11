"""This module stores the CadImporter class, which is used to load CAD parts."""

import math
import os
# needed to prevent dxfgrabber from crashing on import
os.environ['DXFGRABBER_CYTHON'] = 'OFF'
import dxfgrabber # needed for dxf files

from . import geometry
from . import part

class CadImporter(object):
    """Makes an object which can import cad parts.

    Args:
        feamodel (FeaModel): model that we want to import parts into
        layer (int): layer to import, all entities will be flattened
            to one plane. -1 mean all layers, any other number is a specific layer

    Attributes:
        __fea (FeaModel): model
        __fname (str): file we want to import
        __layer (int): layer to import
        __swapxy (bool): If true, swap part from xy to yx orientation
    """

    def __init__(self, feamodel, fname='', layer=-1, swapxy=False):
        self.__fea = feamodel
        self.__fname = fname
        self.__layer = layer
        self.__swapxy = swapxy

    def load(self):
        """Loads the self.__fname cad file

        Returns:
            list: list of Part
        """
        if self.__fname == '':
            print('You must pass in a file name to load!')
            return []
        else:
            fname_list = self.__fname.split('.')
            ext = fname_list[1]
            if ext == 'dxf':
                return self.load_dxf()

    def __fix_point(self, point):
        """Adjusts the point to be in the right plane (yx)"""
        if self.__swapxy:
            return geometry.Point(point.y, point.x)
        return point

    def __get_pt(self, points, point):
        """Returns a point if it is within accuracy of point in points"""
        for realpoint in points:
            dist = point - realpoint
            dist = dist.length()
            if dist < geometry.ACC:
                return realpoint
        return point
    
    def __get_pts_lines(self, lines, arcs):
        """Returns a set of points, and a list of Lines and Arcs

        Args:
            lines (dxfgrabber LINE list): dxf lines
            arcs (dxfgrabber ARC list): dxf arcs

        Returns:
            list: [points, list of Line and Arc]
        """
        # store unique points
        point_set = set()
        all_lines = []
        for ind, line in enumerate(lines):
            start = geometry.Point(line.start[0], line.start[1])
            end = geometry.Point(line.end[0], line.end[1])
            start, end = self.__fix_point(start), self.__fix_point(end)
            point_set.update([start, end])
            line = geometry.Line(start, end)
            all_lines.append(line)
        for ind, arc in enumerate(arcs):
            # dxfgrabber arcs are stored ccw when looking at xy plane
            # x horizontal
            # y vertical
            center = geometry.Point(arc.center[0], arc.center[1])
            sign = 1
            if self.__swapxy == False:
                sign = -1
            center = self.__fix_point(center)
            point_set.add(center)
            startangle = arc.startangle*sign
            endangle = arc.endangle*sign
            angle = endangle - startangle
            if arc.endangle < arc.startangle:
                angle = angle + 360*sign
            """
            print('---------------------------------------')
            print('| ARC')
            print('center: %s' % center)
            print('startangle: %f' % startangle)
            print('endangle: %f' % endangle)
            print('traversed_angle: %f' % angle)
            """
            start_vect = geometry.Point(0, arc.radius)
            if self.__swapxy == False:
                start_vect = geometry.Point(arc.radius, 0)
            start_vect.rot_ccw_deg(arc.startangle*sign)
            end_vect = geometry.Point(0, arc.radius)
            if self.__swapxy == False:
                end_vect = geometry.Point(arc.radius, 0)
            end_vect.rot_ccw_deg(arc.endangle*sign)
            start = center + start_vect
            end = center + end_vect
            start = self.__get_pt(point_set, start)
            end = self.__get_pt(point_set, end)
            point_set.update([start, end])
            rvect = start - center
            if abs(angle) <= 90:
                arc = geometry.Arc(start, end, center)
                all_lines.append(arc)
                #print('1 arc made')
                #print(' %s' % arc)
            else:
                pieces = math.ceil(abs(angle)/90)
                #print('%i arcs being made' % pieces)
                points = [start, end]
                # 2 pieces need 3 points, we have start + end already --> 1 pt
                inserts = pieces + 1 - 2
                piece_ang = angle/pieces
                #print('piece_ang = %f' % piece_ang)
                while inserts > 0:
                    rvect.rot_ccw_deg(piece_ang)
                    point = center + rvect
                    points.insert(-1, point)
                    inserts = inserts - 1
                for ind in range(len(points)-1):
                    point_set.update([points[ind], points[ind+1]])
                    arc = geometry.Arc(points[ind], points[ind+1], center)
                    #print(' %s' % arc)
                    all_lines.append(arc)
        return [list(point_set), all_lines]

    def load_dxf(self):
        """Loads in a dxf file and returns a list of parts

        Returns:
            list: list of Part
        """
        print('Loading file: %s' % self.__fname)
        dwg = dxfgrabber.readfile(self.__fname)
        lines = [item for item in dwg.entities if item.dxftype == 'LINE']
        arcs = [item for item in dwg.entities if item.dxftype == 'ARC']
        if self.__layer > -1:
            lines = [item for item in lines if item.layer == self.__layer]
            arcs = [item for item in arcs if item.layer == self.__layer]
        print('File read.')
        print('Loaded %i lines' % len(lines))
        print('Loaded %i arcs' % len(arcs))

        # get all points and Line and Arc using pycalculix entities
        all_points, all_lines = self.__get_pts_lines(lines, arcs)
        # the index of the point in the set can be used as a hash
        lines_from_ptind = {}
        for line in all_lines:
            ind1 = all_points.index(line.pt(0))
            ind2 = all_points.index(line.pt(1))
            for ind in [ind1, ind2]:
                if ind not in lines_from_ptind:
                    lines_from_ptind[ind] = [line]
                else:
                    if line not in lines_from_ptind[ind]:
                        lines_from_ptind[ind].append(line)

        # make line loops now
        loops = []
        line = all_lines[0]
        this_loop = geometry.LineLoop()
        while len(all_lines) > 0:
            this_loop.append(line)
            if line in all_lines:
                all_lines.remove(line)
            point = line.pt(1)
            ind = all_points.index(point)
            pt_lines = lines_from_ptind[ind]
            if line in pt_lines:
                pt_lines.remove(line)
            if len(pt_lines) == 1:
                # we have the next line
                next_line = pt_lines[0]
                if line.pt(1) != next_line.pt(0):
                    next_line.reverse()
                line = next_line
            elif len(pt_lines) > 1:
                print('One point was connected to > 2 lines.')
                print('Only import simple part loops, or surfaces.')
            if this_loop.closed == True:
                loops.append(this_loop)
                this_loop = geometry.LineLoop()

        # check loops to see if one is inside another
        # if no loops inside self, then self is a part
        # each part in parts is a list of line loops [exterior, hole1, hole2]
        parts = []
        for ind, loop in enumerate(loops):
            other_loops = loops[ind+1:]
            is_exterior = True
            for other_loop in other_loops:
                if loop.inside(other_loop):
                    is_exterior = False
                    break
            if is_exterior:
                # exterior must be clockwise
                if loop.ccw == True:
                    loop.reverse()
                parts.append(loop)
        # remove the exterior loops from our loops list
        for loop in parts:
            loops.remove(loop)
        # now place the child loops under the part exterior loops
        for ind in range(len(parts)):
            exterior_loop = parts[ind]
            holes = []
            # find child holes
            for loop in loops:
                if loop.inside(exterior_loop):
                    # holes must be ccw
                    loop.hole = True
                    if loop.ccw == False:
                        loop.reverse()
                    holes.append(loop)
            loop_list = [exterior_loop] + holes
            parts[ind] = loop_list
            # remove child holes from loop list
            for hole in holes:
                loops.remove(hole)

        # make parts
        parts_list = []
        for part_loops in parts:
            this_part = part.Part(self.__fea)
            for ind, loop in enumerate(part_loops):
                is_hole = loop.hole
                start = loop[0].pt(0)
                this_part.goto(start.x, start.y, is_hole)
                for item in loop:
                    if isinstance(item, geometry.Line):
                        end = item.pt(1)
                        this_part.draw_line_to(end.x, end.y)
                    elif isinstance(item, geometry.Arc):
                        end = item.pt(1)
                        center = item.actr
                        this_part.draw_arc(end.x, end.y, center.x, center.y)
            parts_list.append(this_part)
        print('Parts created: %i' % len(parts_list))
        return parts_list
