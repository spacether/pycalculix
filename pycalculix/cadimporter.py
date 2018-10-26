"""This module stores the CadImporter class, which is used to load CAD parts."""

import collections
import math
import os
import pdb
# needed to prevent dxfgrabber from crashing on import
os.environ['DXFGRABBER_CYTHON'] = 'OFF'
import dxfgrabber # needed for dxf files
import subprocess # needed to run gmsh to make geo files

from . import geometry
from . import partmodule

class CadImporter(object):
    """Makes an object which can import cad parts.

    Args:
        feamodel (FeaModel): model that we want to import parts into
        layer (int): layer to import, all entities will be flattened
            to one plane. -1 mean all layers, any other number is a specific layer
        swapxy (bool): True rotates the part from x axial to y axial
        scale (str): string telling the unit conversion scalar

            * Any option of 'fromunit-tounit' using the below units
            * mm, m, in, ft
            * Examples 'mm-m' 'm-in' 'ft-mm' 'in-ft'
            * Default value is '' and does not apply a scale factor

    Attributes:
        __fea (FeaModel): model
        __fname (str): file we want to import
        __layer (int): layer to import
        __swapxy (bool): If true, swap part from xy to yx orientation
    """

    def __init__(self, feamodel, fname='', layer=-1, swapxy=False, scale=''):
        self.__fea = feamodel
        self.__fname = fname
        self.__layer = layer
        self.__swapxy = swapxy
        self.__scale = scale

    def load(self):
        """Loads the self.__fname cad file

        Returns:
            list: list of Part
        """
        if self.__fname == '':
            print('You must pass in a file name to load!')
            return []

        ext = os.path.splitext(self.__fname)[1]
        first_pt = None
        if len(self.__fea.points) > 0:
            first_pt = self.__fea.points[0]
        if ext == '.dxf':
            parts = self.__load_dxf()
        elif ext in ['.brep', '.brp', '.iges', '.igs', '.step', '.stp']:
            self.__make_geo()
            parts = self.__load_geo()
        last_pt = None
        if first_pt != None:
            if len(self.__fea.points) > 2:
                last_pt = self.__fea.points[-1]
        if self.__scale != '':
            # call scale
            pass
        return parts

    def __fix_tuple(self, xy_tup):
        """Adjusts the point to be in the right plane (yx)"""
        if self.__swapxy:
            return xy_tup[::-1]
        return xy_tup

    @staticmethod
    def __find_make_pt(xy_tup, points_dict):
        """
        Returns a point if it exists within geometry.ACC
        or makes, stores and returns the point if it doesn't exist
        """
        point = points_dict.get(xy_tup)
        if point is not None:
            return point
        xy_point = geometry.Point(xy_tup[0], xy_tup[1])
        for realpoint in points_dict.values():
            dist = (xy_point - realpoint).length()
            if dist < geometry.ACC:
                return realpoint
        points_dict[xy_tup] = xy_point
        return xy_point

    def __get_pts_lines(self, lines, arcs):
        """Returns a set of points, and a list of Lines and Arcs

        Args:
            lines (dxfgrabber LINE list): dxf lines
            arcs (dxfgrabber ARC list): dxf arcs

        Returns:
            list: [list of points, list of Line and Arc]
        """
        # store unique points
        points_dict = {}
        all_lines = []
        for ind, line in enumerate(lines):
            tup = self.__fix_tuple((line.start[0], line.start[1]))
            start = self.__find_make_pt(tup, points_dict)
            tup = self.__fix_tuple((line.end[0], line.end[1]))
            end = self.__find_make_pt(tup, points_dict)
            line = geometry.Line(start, end)
            all_lines.append(line)
        for ind, arc in enumerate(arcs):
            # dxfgrabber arcs are stored ccw when looking at xy plane
            # x horizontal
            # y vertical
            tup = self.__fix_tuple((arc.center[0], arc.center[1]))
            center = self.__find_make_pt(tup, points_dict)
            sign = -1
            if self.__swapxy:
                sign = 1
            startangle = arc.start_angle*sign
            endangle = arc.end_angle*sign
            angle = endangle - startangle
            if arc.end_angle < arc.start_angle:
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
            start_vect.rot_ccw_deg(arc.start_angle*sign)
            end_vect = geometry.Point(0, arc.radius)
            if self.__swapxy == False:
                end_vect = geometry.Point(arc.radius, 0)
            end_vect.rot_ccw_deg(arc.end_angle*sign)
            start = center + start_vect
            start_tup = (start.x, start.y)
            end = center + end_vect
            end_tup = (end.x, end.y)
            start = self.__find_make_pt(start_tup, points_dict)
            end = self.__find_make_pt(end_tup, points_dict)
            rvect = start - center
            if abs(angle) <= 90:
                arc = geometry.Arc(start, end, center)
                all_lines.append(arc)
                print('1 arc made')
                continue
                #print(' %s' % arc)
            pieces = math.ceil(abs(angle)/90)
            print('%i arcs being made' % pieces)
            points = [start, end]
            # 2 pieces need 3 points, we have start + end already --> 1 pt
            inserts = pieces + 1 - 2
            piece_ang = angle/pieces
            #print('piece_ang = %f' % piece_ang)
            while inserts > 0:
                rvect.rot_ccw_deg(piece_ang)
                point = center + rvect
                tup = (point.x, point.y)
                point = self.__find_make_pt(tup, points_dict)
                points.insert(-1, point)
                inserts = inserts - 1
            for ind in range(len(points)-1):
                #print(' %s' % arc)
                arc = geometry.Arc(points[ind], points[ind+1], center)
                all_lines.append(arc)
        for line in all_lines:
            line.save_to_points()
        return [list(points_dict.values()), all_lines]

    def __make_geo(self):
        """Makes a gmsh geo file given a step, iges, or brep input"""
        # gmsh freecad_part.iges -o out_iges.geo -0
        fname_list = self.__fname.split('.')
        geo_file = fname_list[0]+'.geo'
        runstr = "%s %s -o %s -0" % (environment.GMSH, self.__fname, geo_file)
        print(runstr)
        subprocess.call(runstr, shell=True)
        print('Wrote file: %s' % geo_file)

    def __load_geo(self):
        """Loads in a gmsh geo file and returns a list of parts

        Returns:
            list: list of Part
        """
        pass
        # process any splines? and turn them into arcs
        # http://www.mathopenref.com/constcirclecenter.html
        # find max dist between points
        # double it
        # select two segments
        # draw normal lines
        # find intersections, that is the center

    @staticmethod
    def __dangling_points(all_points):
        return [point for point in all_points
                if len(point.lines) == 1 and not point.arc_center]

    def __load_dxf(self):
        """Loads in a dxf file and returns a list of parts

        Returns:
            list: list of Part
        """
        # pdb.set_trace()
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
        print('Loaded %i line segments, lines or arcs' %
              (len(lines)+len(arcs)))
        # get all points and Line and Arc using pycalculix entities
        print('Converting to pycalculix lines arcs and points ...')
        all_points, all_lines = self.__get_pts_lines(lines, arcs)
        print('Loaded %i line segments, lines or arcs' % len(all_lines))
        print('Loaded %i points' % len(all_points))
        # for point in all_points:
        #     print('%s %s' % (point, point.lines))
        # for line in all_lines:
        #     print('%s %s' % (line, line.points))

        # remove all lines that are not part of areas
        dangling_points = self.__dangling_points(all_points)
        # pdb.set_trace()
        pruned_geometry = bool(dangling_points)
        while dangling_points:
            for point in dangling_points:
                all_points.remove(point)
                print('Removed point= %s' % point)
                dangling_line = list(point.lines)[0]
                point.unset_line(dangling_line)
                if dangling_line in all_lines:
                    all_lines.remove(dangling_line)
                    print('Removed line= %s' % dangling_line)
            dangling_points = self.__dangling_points(all_points)
        if pruned_geometry:
            print('Remaining line segments: %i' % len(all_lines))
            print('Remaining points: %i' % len(all_points))

        # make line all_loops now
        all_loops = []
        line = all_lines[0]
        this_loop = geometry.LineLoop()
        while len(all_lines) > 0:
            this_loop.append(line)
            all_lines.remove(line)
            if this_loop.closed == True:
                all_loops.append(this_loop)
                this_loop = geometry.LineLoop()
                if all_lines:
                    line = all_lines[0]
                continue
            point = line.pt(1)
            other_lines = point.lines - set([line])
            if len(other_lines) > 1:
                # note: one could exclude connected segment nodes
                # make disconnected line all_loops, then have another
                # loop to connect thos disconnected line all_loops
                print('One point was connected to > 2 lines.')
                print('Only import simple part all_loops, or surfaces.')
                raise Exception('Import geometry is too complex')
            next_line = list(other_lines)[0]
            if line.pt(1) != next_line.pt(0):
                next_line.reverse()
            line = next_line

        # find exterior loops
        exterior_loops = []
        for ind, loop in enumerate(all_loops):
            other_loops = all_loops[ind+1:]
            other_loops.extend(exterior_loops)
            is_exterior = True
            for other_loop in other_loops:
                if loop.inside(other_loop):
                    is_exterior = False
                    break
            if is_exterior:
                # exterior must be clockwise
                if loop.ccw:
                    loop.reverse()
                exterior_loops.append(loop)
        # remove the found part exterior loops from all_loops
        for exterior_loop in exterior_loops:
            all_loops.remove(exterior_loop)
        # each part in parts is a list of line all_loops
        # [exterior, hole1, hole2]
        parts = [[exterior_loop] for exterior_loop in exterior_loops]
        # now place the child hole loops after the part exterior loop
        for part_loops in parts:
            exterior_loop = part_loops[0]
            # find child holes
            for hole_loop in all_loops:
                if hole_loop.inside(exterior_loop):
                    hole_loop.hole = True
                    # holes must be ccw
                    if not hole_loop.ccw:
                        hole_loop.reverse()
                    part_loops.append(hole_loop)
            # remove child holes from loop list
            for hole_loop in part_loops[1:]:
                all_loops.remove(hole_loop)

        # make parts
        parts_list = []
        for part_loops in parts:
            this_part = partmodule.Part(self.__fea)
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
