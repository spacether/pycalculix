#!/usr/bin/python
"""
calculix.py
Description:
This module provides a python interface to the calculix FEA preprocessor and
solver. Simplicity is emphasized.
Axisymmetric, plane stress and plane strain elements are supported.
First and second order triangles and quadrilaterals are supported.
The python file stores set-up information like geometry and loads, then
produces the required cgx mesh building file (fbd), and uses the resultant
cgx mesh and load files to produce the ccx input (inp) file.
The cgx preprocessor and the ccx solver are run in the background,
allowing for fast iteration and answers.

Usefull applications:
-Quick Kt analysis of simple 2D geometry.
-Trade studies for flat plate parts or axisymmetric parts

Notes:
I build a chunker in python which tries to cut big areas (> 5 sides) which
cgx can't mesh into smaller areas (<= 5 sides) which are meshable in cgx.
The chunker may not always work. I'd like to integrate gmsh in the future
so mesh quality can be more controllable, and less error prone.

Future Goals:
-Gmsh integration including holes in parts
-Contact support between parts
  -2d axisymmetric bolted joint modeling
  -2d axisymmetric interferences (axial or radial)
-Openfoam solver runs
-Mapping open foam results back on to calculix models
-Rocket nozzle modeling, struct, thermal, and cfd
"""

__author__ = "Justin Black"
__copyright__ = "Copyright 2014, Justin Black"
__credits__ = ["Justin Black"]
__license__ = "GPL V2"
__version__ = "0.9"
__maintainer__ = "Justin Black"
__email__ = "justin.a.black[incert~at~sine]gmail.com"
__status = "Production"

# this makes, nodes, lines, and an area(s), and part(s)
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon  # needed for plotting elements
from matplotlib.collections import PatchCollection  # element plotting
from matplotlib.patches import Arc as AArc
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.colorbar as colorbar      

from math import atan2, pi, cos, sin, radians, ceil, acos, sqrt
import os # used to delete files written by cgx
import subprocess # used to launch cgx and ccx
import re # used to get info from frd file
from numpy import roots # needed to find principal stresses
import platform # need this to check for windows and do the below dpi fix for high dpi monitors

# set the dpi if the user is running Windows 8 or higher
_dpi = None
if platform.system() == 'Windows':
    if float( platform.release() ) >= 8:
        ''' only run if we're on windows 8 or higher
        I need to reset the default DPI for my high DPI monitor because
        matplotlib does not set the DPI higher when it should
        '''
    
        # Get windows dpi:
        import ctypes
        user32 = ctypes.windll.user32
        w_curr = user32.GetSystemMetrics(0)
        user32.SetProcessDPIAware()
        w_phys = user32.GetSystemMetrics(0)
        _dpi = round(w_phys*96/w_curr,0)
        
        from pylab import rcParams
        rcParams['figure.dpi'] = _dpi

#accuracy for small numbers in math below
_acc = .00001


_ccx = r'C:\Program Files (x86)\bConverged\CalculiX\ccx\ccx.exe'
_cgx = r'C:\Program Files (x86)\bConverged\CalculiX\cgx\cgx.exe'
_gmsh = r'C:\Program Files (x86)\gmsh-2.8.5\gmsh.exe'


# Calculix CGX elements
# axisymmetric
_cgx_elements = {}
_cgx_elements['tri2axisym'] = 'tr6c'
_cgx_elements['tri1axisym'] = 'tr3c'
_cgx_elements['quad2axisym'] = 'qu8c'
_cgx_elements['quad1axisym'] = 'qu4c'
# plane stress
_cgx_elements['tri2plstress'] = 'tr6s'
_cgx_elements['tri1plstress'] = 'tr3s'
_cgx_elements['quad2plstress'] = 'qu8s'
_cgx_elements['quad1plstress'] = 'qu4s'
# plane strain
_cgx_elements['tri2plstrain'] = 'tr6e'
_cgx_elements['tri1plstrain'] = 'tr3e'
_cgx_elements['quad2plstrain'] = 'qu8e'
_cgx_elements['quad1plstrain'] = 'qu4e'

# Calculix CCX elements
# axisymmetric
_ccx_elements = {}
_ccx_elements['tri2axisym'] = 'CAX6'
_ccx_elements['tri1axisym'] = 'CAX3'
_ccx_elements['quad2axisym'] = 'CAX8'
_ccx_elements['quad1axisym'] = 'CAX4'
# plane stress
_ccx_elements['tri2plstress'] = 'CPS6'
_ccx_elements['tri1plstress'] = 'CPS3'
_ccx_elements['quad2plstress'] = 'CPS8'
_ccx_elements['quad1plstress'] = 'CPS4'
# plane strain
_ccx_elements['tri2plstrain'] = 'CPE6'
_ccx_elements['tri1plstrain'] = 'CPE3'
_ccx_elements['quad2plstrain'] = 'CPE8'
_ccx_elements['quad1plstrain'] = 'CPE4'

# element colors, 0-1, 0=black, 1=whate
_ecolor = '.4'
_fcolor = '.9'


# constants for results files
_resfields = {}
_resfields['displ'] = 'ux,uy,uz,utot'.split(',')
_resfields['stress'] = 'Sx,Sy,Sz,Sxy,Syz,Szx,Seqv,S1,S2,S3'.split(',')
_resfields['strain'] = 'ex,ey,ez,exy,eyz,ezx,eeqv,e1,e2,e3'.split(',')
_resfields['force'] = 'fx,fy,fz'.split(',')
# make an inverse dict where passing 'ux' gives us 'displ'
_fieldtype = {}
for (k, v) in _resfields.items():
    for vi in v:
        _fieldtype[vi] = k

def frange(start, end=None, inc=None):
    """A range function, that does accept float increments..."""

    if end == None:
        end = start + 0.0
        start = 0.0
    else: start += 0.0 # force it to be a float

    if inc == None:
        inc = 1.0
    else:
        if end < start and inc > 0:
            inc = inc*-1.0
    count = int(ceil((end - start) / inc)) + 1

    L = [None,] * count

    L[0] = start
    for i in range(1,count):
        L[i] = L[i-1] + inc
    return L

class Idobj():
    # this stores an id number for nodes, lines, areas etc
    def __init__(s):
        s.id = -1
    def set_id(s, id):
        s.id = id    
    
class Point(Idobj):
    # this stores a point or vector
    def __init__(s, x, y):
        s.x = x
        s.y = y
        Idobj.__init__(s)
    def get_name(s):
        # returns item name
        return 'P'+str(s.id)
    def __eq__(s, other):
        if s.x == other.x and s.y == other.y:
            return True
        return False
    def __add__(s, other):
        # this doesn't need to return point with an ID
        return Point(s.x+other.x,s.y+other.y)
    def __sub__(s, other):
        # this doesn't need to return point with an ID
        return Point(s.x-other.x,s.y-other.y)
    def __mul__(s, factor):
        # this doesn't need to return point with an ID
        return Point(s.x*factor, s.y*factor)
    def length(s):
        # returns distance from this point to 0,0
        res = (s.x**2 + s.y**2)**0.5
        return res
    def make_unit(s):
        # converts the current vector into a unit vector
        lenval  = s.length()
        s.x = s.x*1/lenval
        s.y = s.y*1/lenval
    def ang_rad(s):
        # returns angle, assume point is vector from 0,0
        vert = s.x
        horiz = s.y
        radians = atan2(vert, horiz)
        return radians
    def ang_deg(s):
        # returns angle, assume point is vector from 0,0
        radians = s.ang_rad()
        deg = radians * 180.0 / pi
        return deg
    def rot_ccw_deg(s, ang):
        # this rotates the current vector by ccw degrees about 0,0
        ang = radians(ang)
        ax = s.y*cos(ang) - s.x*sin(ang)
        rad = s.y*sin(ang) + s.x*cos(ang)
        ax = round(ax, 5)
        rad = round(rad, 5)
        (s.x, s.y) = (rad, ax)
    def ang_bet_rad(s, other):
        # returns the angle between two vectors
        # self is the first vector
        # other is the second
        avect = Point(s.y*other.x - s.x*other.y, s.y*other.y + s.x*other.x)
        ang = avect.ang_rad()
        return ang
    def ang_bet_deg(s, other):
        # returns the angle between two vectors
        # self is the first vector
        # other is the second
        ang1 = s.ang_deg()
        ang2 = other.ang_deg()
        avect = Point(s.y*other.x - s.x*other.y, s.y*other.y + s.x*other.x)
        ang = avect.ang_deg()
        #print('a1=%f, a2=%f, delta=%f' % (ang1, ang2, ang))
        return ang
    def __str__(s):
        val = 'Point, id %i, (x,y)=(%f,%f)' % (s.id, s.x,s.y)
        return val


class Line(Idobj):
    def __init__(s, p1, p2, sign=1):
        s._points_ = [p1, p2]
        s.sign = sign
        s.midpt = s.mid()
        Idobj.__init__(s)
    def get_name(s):
        # returns line name
        return 'L'+str(s.id)
    def plot(s, ax):
        # plots the line on passes matplotlib ax
        lax = [pt.y for pt in s._points_]
        lrad = [pt.x for pt in s._points_]
        ax.plot(lax, lrad)
        pax = s.midpt.y
        prad = s.midpt.x
        ax.annotate(s.get_name(), (pax,prad))
    def length(s):
        # returns the length of the line
        return (s.pt(0)-s.pt(1)).length()
    def set_ediv(s, ediv):
        # this sets the number of dividions on the line
        s.ediv = ediv
    def reverse(s):
        # makes a copy of this line but fips the sign
        a = Line(s._points_[0], s._points_[1], s.sign*-1)
        a.id = s.id
        return a
    def pt(s, ind):
        # returns points in correct order, pass it 0 or 1
        if s.sign == -1:
            return s._points_[not ind]
        else:
            return s._points_[ind]        
    def set_pt(s, ind, pt):
        # used to update specific points in the line
        if s.sign == -1:
            s._points_[not ind] = pt
        else:
            s._points_[ind] = pt        
        s.midpt = s.mid()
    def mid(s):
        pt = (s.pt(0) + s.pt(1))*0.5
        return pt
    def touches(s, other):
        # this checks if a line is touching the passed line
        [p1, p2] = [s.pt(0), s.pt(1)]
        if p1 in other._points_ or p2 in other._points_:
            return True
        else:
            return False
    def offset(s, dist):
        # this returns a temperary line which is offset from the line
        lvect = s.pt(1) - s.pt(0)
        # rotate it ccw by 90 degrees, and make it a unit
        lvect.rot_ccw_deg(90)
        lvect.make_unit()
        # scale it by the new amount
        lvect = lvect*dist
        # add the vector onto the defining points
        p0 = s.pt(0) + lvect
        p1 = s.pt(1) + lvect
        tmpline = Line(p0, p1)
        return tmpline
    def get_abc(s):
        # this returns a list of abc terms for a line definition
        # ax + by + c = 0
        lpt = s.pt(1)-s.pt(0)
        dx = lpt.x 
        dy = lpt.y
        (a, b, c) = (0.0, 0.0, 0.0)
        if dy == 0:
            # vertical radial line
            (b, c) = (1.0, s.pt(0).y*-1.0)
        elif dx == 0:
            # horizontal axial line
            (a, c) = (1.0, s.pt(0).x*-1.0)
        else:
            slope = dy/dx
            offset = s.pt(0).y - slope*s.pt(0).x
            (a, b, c) = (1.0, -1.0*slope, -1.0*offset)
        return [a,b,c]
    def get_perp_vec(s, pt=None):
        # returns vector perpendicular to current one
        # pt is not needed for line but is needed for arc
        lvect = s.pt(1) - s.pt(0)
        # rotate it ccw by 90 degrees
        lvect.rot_ccw_deg(90)
        return lvect
    def get_abc_perp(s, pt):
        # this returns a list of abc terms for a line definition
        # ax + by + c = 0
        # make a vector of the current line
        lvect = s.get_perp_vec()
        # make the perpendicular line end point
        p1 = pt + lvect
        tmpline = Line(pt, p1)
        return tmpline.get_abc()
    def arc_tang_intersection(s, pt, mag):
        # returns a point of intersection on this line
        # by drawing a perp line passing through arc center pt
        # this is used to trim lines when we put a fillet in

        v = s.get_perp_vec()
        v.make_unit()
        p2 = pt+v*(-2*mag)

        tmpline = Line(pt, p2)
        newpt = tmpline.intersects(s)
        if type(newpt) == type(None):
            print('Intersection failed!')
        return newpt

    def coincident(s, pt):
        # checks if the given point is on this line
        # ax + by + c = 0
        [a,b,c] = s.get_abc()
        remainder = a*pt.x + b*pt.y + c
        if remainder < _acc:
            # point is on the line equation, but is it between end pts
            lvect = s.pt(1) - s.pt(0)
            # point0 + lvect*term = pt
            nondim = pt - s.pt(0)
            nondim = nondim.x/lvect.x
            if 0 <= nondim <= 1.0:
                # we're on the line
                return True
            else:
                # we're off of the line
                return False
        else:
            # point is not on line equation
            return False
    def intersects(s, other):
        # returns intersection point between this line and other line
        # other is a straight line
        if isinstance(other, Line):
            slope1 = s.get_abc()
            slope2 = other.get_abc()
            if slope1[:2] == slope2[:2]:
                # lines are paralel and will not have an intersection point
                return None
            else:
                # parametric intersection calculation
                # p1 -> p2 = ta parametric, p2 -> p3 tb parametric
                (p1, p2, p3, p4) = (s.pt(0), s.pt(1), other.pt(0), other.pt(1))
                (x1,x2,x3,x4)= (p1.y, p2.y, p3.y, p4.y)
                (y1,y2,y3,y4)= (p1.x, p2.x, p3.x, p4.x)
                (x31,x21,x43)=(x3-x1,x2-x1,x4-x3)
                (y31,y21,y43)=(y3-y1,y2-y1,y4-y3)
                tanum = x43*y31 - x31*y43
                tadenom = x43*y21 - x21*y43
                tbnum = x21*y31 - x31*y21
                tbdenom = x43*y21 - x21*y43
                if tadenom != 0 and tbdenom != 0:
                    ta = tanum/tadenom
                    tb = tbnum/tbdenom
                    if 0 <= ta <= 1.0 and 0 <= tb <= 1.0:
                        pnew = (p2-p1)
                        pt = p1 + pnew*ta
                        return pt
                    else:
                        # intersection is outside the line bounds
                        return None
                else:
                    return None
        else:
            # arc line intersection
            return None
    def __str__(s):
        p0 = s.pt(0)
        p1 = s.pt(1)
        val = 'Line, id %i p0: %s p1: %s' % (s.id, p0, p1)
        return val        

class Arc(Idobj):
    # this arc should be 90 degrees max
    def __init__(s, p1, p2, actr, sign=1):
        s._points_ = [p1, p2]
        s.sign = sign
        s.actr = actr
        s.radius = (p1-actr).length()
        s.concavity = s.get_concavity()
        s.midpt = s.mid()
        Idobj.__init__(s)        
    def get_name(s):
        # returns line name
        return 'L'+str(s.id)
    
    def plot(s, ax):
        # plot the arc on passed ax
        (pax, prad) = (s.midpt.y, s.midpt.x)
        ax.annotate(s.get_name(), (pax,prad))
        
        # need center, radius, ang1, and ang2
        ctr = s.actr
        vect1 = None
        sign = 1
        if s.concavity == 'concave':
            vect1 = s.pt(0)-ctr
        else:
            sign = -1
            vect1 = s.pt(1)-ctr
        rad = s.radius
        ang1 = (vect1.ang_deg())
        ang2 = ang1 + sign*s.get_ang()

        # matplotlib assumes ccw arc drawing, calculix assumes cw drawing, reverse order
        a = AArc(xy=[ctr.y, ctr.x], width=2*rad, height=2*rad, angle=0, theta1=ang1, theta2=ang2)
        ax.add_artist(a)

    def length(s):
        # returns the length of the arc
        a = s.pt(0)-s.actr
        b = s.pt(1)-s.actr
        ang = a.ang_bet_rad(b)
        res = s.radius*ang
        return res
    def set_ediv(s, ediv):
        # this sets the number of dividions on the line
        s.ediv = ediv
    def reverse(s):
        # makes a copy of this line but fips the sign
        a = Arc(s._points_[0], s._points_[1], s.sign*-1)
        a.id = s.id
        return a
    def pt(s, ind):
        # returns points in correct order
        if s.sign == -1:
            return s._points_[not ind]
        else:
            return s._points_[ind]        
    def mid(s):
        # this returns the midpoint of the arc 
        pt = s.get_pt_at(0.5)
        return pt
    def touches(s, other):
        # checks if this arc touches the passed other line
        [p1, p2] = [s.pt(0), s.pt(1)]
        if p1 in other._points_ or p2 in other._points_:
            return True
        else:
            return False
    def sgn(s, val):
        # helper function for interection code
        if val < 0:
            return -1
        else:
            return 1
    def get_ang(s):
        # returns angle or ray from beginning to end of arc from arc ctr in deg

        # answer is between 180 and -180
        # 40 to -40, --> -80 (neg is cw), convex r
        # -40 to 40, --> 80, concave l
        # -160 to 160, --> -40 convex l
        # 160 to -160, --> 40, concave r
        # pos is concave
        # neg is convex

        a = s.pt(0)-s.actr
        b = s.pt(1)-s.actr
        ang = a.ang_bet_deg(b)
        return ang
    def get_concavity(s):
        # returns concave or convex
        ang = s.get_ang()        
        res = 'concave'
        if ang < 0:
            res = 'convex'
        return res
        
    def get_pt_at(s, nondim):
        # nondim is a number between 0 and 1, a pt on the arc is returned
        a = s.pt(0)-s.actr

        # neg is convex, pos is concave
        ang = s.get_ang()
        
        ang_delta = nondim*ang
        a.rot_ccw_deg(ang_delta)

        res = s.actr + a
        return res
    def get_perp_vec(s, pt):
        # returns perp vector
        resv = None
        if s.concavity == 'convex':
            resv = pt - s.actr
        elif s.concavity == 'concave':
            resv = s.actr - pt
        return resv
    def get_abc_perp(s, pt):
        # this returns a list of abc terms for a line perp to the arc on pt
        # ax + by + c = 0
        
        # if the arc is convex, use the vector pt-actr
        # if the arc is concave, use the actr - pt        
        resv = s.get_perp_vec(pt)
        (dx, dy) = (resv.x, resv.y)        

        (a, b, c) = (0.0, 0.0, 0.0)
        if dy == 0:
            # vertical radial line
            (b, c) = (1.0, pt.y*-1.0)
        elif dx == 0:
            # horizontal axial line
            (a, c) = (1.0, pt.x*-1.0)
        else:
            slope = dy/dx
            offset = pt.y - slope*pt.x
            (a, b, c) = (1.0, -1.0*slope, -1.0*offset)
        return [a,b,c]

    def coincident(s, pt):
        # checks if the given point is on this arc
        # (x - xc)^2 + (y - yc)^2 - r^2 = 0
        remainder = (pt.x - s.actr.x)**2 + (pt.y - s.actr.y)**2 - s.radius**2
        if remainder < _acc:
            # point is on the circle, but is it in the arc?
            # check if the angle traversed to the point is within the ang
            # traversed by the arc
            a1 = s.pt(0)-s.actr
            a2 = pt-s.actr
            ang_pt = a1.ang_bet_deg(a2)
            
            # this is the angle traversed byt the arc
            ang = s.get_ang()
            minval = min(0, ang)
            maxval = max(0, ang)
            if minval <= ang_pt <= maxval:
                # point is on our arc
                return True
            else:
                # point on circle but outside our arc
                return False
        else:
            # point not on arc
            return False

    def intersects(s, other):
        if isinstance(other, Line):
            # arc-line intersection
            # reset the math to be about the circle centroid
            # formula from: http://mathworld.wolfram.com/Circle-LineIntersection.html
            p1 = other.pt(0) - s.actr
            p2 = other.pt(1) - s.actr
            # for us x and y are reversed because of y axial
            (x1, y1, x2, y2) = (p1.y, p1.x, p2.y, p2.x)
            dx = x2 - x1
            dy = y2 - y1
            r = s.radius
            dr = (dx**2 + dy**2)**(0.5)
            D = x1*y2 - x2*y1
            discr = 1.0*(r**2)*(dr**2) - D**2
            if discr < 0:
                # no intersection
                return None
            else:
                # tangent or intersection
                x1 = (D*dy+s.sgn(dy)*dx*((r**2)*(dr**2)-D**2)**(0.5))/(dr**2)
                y1 = (-D*dx+abs(dy)*dx*((r**2)*(dr**2)-D**2)**(0.5))/(dr**2)
                x2 = (D*dy-s.sgn(dy)*dx*((r**2)*(dr**2)-D**2)**(0.5))/(dr**2)
                y2 = (-D*dx-abs(dy)*dx*((r**2)*(dr**2)-D**2)**(0.5))/(dr**2)

                # convert result points back to global coordinates
                p1 = Point(y1, x1) + s.actr
                p2 = Point(y2, x2) + s.actr
                res = []
                # check that the resultant points are on line and arcs
                if s.coincident(p1) and other.coincident(p1):
                    res.append(p1)
                if p1 != p2 and s.coincident(p2) and other.coincident(p2):
                    res.append(p2)
                if len(res)==1:
                    return res[0]
                else:
                    return res
        elif isinstance(other, Arc):
            # arc-arc intersection
            # I have not yet written the code here
            pass

    def __str__(s):
        val = 'Arc, id %i' % (s.id,)
        return val

class Area(Idobj):
    def __init__(s, p, line_list=[]):
        # initialtes the area
        s.p = p # parent part reference
        s.closed = False
        s.update(line_list)
        Idobj.__init__(s)
        s.matl = None
        s.etype = None

    def set_etype(s, etype):
        s.etype = etype

    def get_name(s):
        # returns area name
        return 'A'+str(s.id)

    def update(s, line_list):
        # adds the line list and checks for closure
        s.lines = line_list
        if len(line_list) > 2:
            if s.lines[0].pt(0) == s.lines[-1].pt(1):
                s.closed = True
                s.calc_center()        

    def line_from_startpt(s, pt):
        # returns the line that starts on the given point
        for line in s.lines:
            if line.pt(0) == pt:
                return line
        return None

    def get_points(s):
        # gets all the points in the area, does not dupllicate the last point
        pts = []
        for line in s.lines:
            p1 = line.pt(0)
            p2 = line.pt(1)
            if p1 not in pts:
                pts.append(p1)
            if p2 not in pts:
                pts.append(p2)
        return pts

    def calc_center(s):
        # this function calculates the area, and the area centroid
        
        # get list of points, append first point as last
        pts = s.get_points()
        pts.append(pts[0])
        
        # loop to get area and centroid
        cx = 0.0
        cy = 0.0
        area = 0.0
        for i in range(len(pts)-1):
            term1 = (pts[i].x*pts[i+1].y)-(pts[i+1].x*pts[i].y)
            area += 1.0*term1
            cx += (pts[i].x + pts[i+1].x)*term1
            cy += (pts[i].y + pts[i+1].y)*term1
        area = 0.5*area
        cx = (1/(6.0*area))*cx
        cy = (1/(6.0*area))*cy
        s.area = area
        s.centroid = Point(cx, cy)

    def get_maxlength(s):
        # returns the max length between points in an area
        pts = s.get_points()
        maxlen = 0.0
        # loop through points checking dist to next point
        for ind, p1 in enumerate(pts[:-1]):
            for p2 in pts[ind:]:
                vect = p1 - p2
                dist = vect.length()
                if dist > maxlen:
                    maxlen = dist
        return maxlen

    def add(s, line):
        # this adds a line to the area
        s.lines.append(line)
        # check to see if the area is closed
        if line.pt(1) == s.lines[0].pt(0):
            s.closed = True
            s.calc_center()

    def pt_inside(s, pt):
        # returns true if pt is inside area
        # check bounding box first
        # we assume that the point is not on the boundary lines
        pts = s.get_points()
        xs = [p.x for p in pts]
        ys = [p.y for p in pts]
        if (min(xs) <= pt.x <= max(xs)) and (min(ys) <= pt.y <= max(ys)):
            # point is in area bounding box
            vlen = s.get_maxlength()
            pend = Point(pt.x + 0.0001*vlen, pt.y+vlen)
            cutline = Line(pt, pend)
            overlaps = 0
            for line in s.lines:
                ipts = line.intersects(cutline)
                if type(ipts) == type([]) or isinstance(ipts, Point):
                    if type(ipts) == type([]):
                        # we intersected an arc and hit it twice
                        overlaps += 2
                    else:
                        overlaps += 1
            # if overlaps is odd, we are in the polygon
            if overlaps % 2 == 1:
                return True
            else:
                return False
        else:
            # point is outside of area bounding box
            return False
        
    def line_insert(s, lgiven, lnew):
        # this inserts a line after the given line
        if lgiven in s.lines:
            ind = s.lines.index(lgiven)+1
            s.lines.insert(ind,lnew)
            return True
        else:
            return False

# this class is used for lists of nodes, lines, areas, parts, etc
class Item_List(list):
    def __init__(s):
        super().__init__() # these lists start empty

    def get_ids(s):
        # returns list of ids
        return [a.id for a in s]

    def get_next_id(s):
        # returns the next id to use
        ids = s.get_ids()
        minid = 0
        if len(ids) == 0:
            return minid    # list is empty so return the minid number
        else:
            ids = sorted(ids)
            maxid = ids[-1]
            unused = list(set(list(range(minid, maxid+2))) - set(ids))
            return unused[0]    # return first available id number

    def append(s, item):
        # add item to the list and set it's id
        idnum = s.get_next_id()
        item.set_id(idnum)
        super().append(item)
        return item

class ID_List(list):
    def __init__(s):
        list.__init__([])
    def get_minid(s):
        # check for min id
        ids = [x.id for x in s]
        return min(ids)
    def get_maxid(s):
        # returns the max id any item in the list
        ids = [x.id for x in s]
        return max(ids)
    def idget(s, idnum):
        # returns item given an id
        for item in s:
            if item.id == idnum:
                return item
        #print ('ID %i was not found!' % (idnum))
        return None
    def set_minid(s, val):
        print('Re-indexing elements to min element number of: %i' % val)
        print('Old min:%i  max:%i' % (s.get_minid(),s.get_maxid()))
        minid = s.get_minid()
        offset = minid - val
        for item in s:
            item.id -= offset
        print('New min:%i  max:%i' % (s.get_minid(),s.get_maxid()))

def chunk_list(inlist, size):
    # returns a list of lists size <= size
    res = []
    numlists = ceil(len(inlist)/size)
    for ind in range(numlists):
        res.append( inlist[ind*size:(ind+1)*size] )
    return res

class Component(Idobj):
    def __init__(s, item_list, ctype, cname=''):
        # stores components for application of pressure, displacement, etc
        # item list can be a list of points or lines
        # ctype can be 'nall'=all nodes, 'n1'=element 1st order nodes
        #       or 'faces' or 'elements' for elements
        s.items = item_list
        s.ctype = ctype
        s.name = cname+'_'+ctype
        Idobj.__init__(s)
    def set_name(s, name):
        s.name = name
    def get_name(s):
        # returns comp name by id only
        return 'C'+str(s.id)
    def get_children(s):
        # returns nodes or faces under items
        res = []
        for L in s.items:
            children = getattr(L, s.ctype)
            if isinstance(children, list):
                res += children
            else:
                res += [children]
        return res
    def ccx(s):
        # we can only write components of nodes and elements
        res = []
        items_per_line = 6
        firstline = ''
        if s.ctype in ['nodes', 'nall']:
            # node componenet
            firstline = '*NSET,NSET='+s.name
        elif s.ctype == 'elements':
            # element componet
            firstline = '*ELSET,ELSET='+s.name
        res.append(firstline)
        items = s.get_children()
        grouped_items = chunk_list(items, items_per_line)
        for group in grouped_items:
            item_ids = [str(x.id) for x in group]
            line = ', '.join(item_ids)
            if group != grouped_items[-1]:
                line += ','
            res.append(line)
        return res
    
    def write_cgx(s):
        # this returns a list of lines that define and write the component in
        # the preprocessor
        res = []
        
        # set the type of component
        parent = 'l'
        if isinstance(s.items[0], Line) or isinstance(s.items[0], Arc):
            parent = 'l'
        elif isinstance(s.items[0], Point):
            parent = 'p'
        
        # pull the list of entities
        alist = ' '.join([a.get_name() for a in s.items])
        
        if s.ctype == 'n1':
            # write component of line element first order nodes
            res.append('seta '+s.name+' '+parent+' '+alist)                       
        elif s.ctype == 'nall':
            res.append('seta '+s.name+' '+parent+' '+alist)                       
            res.append('comp '+s.name+' do')                       
            # write component of all line element nodes
        elif s.ctype == 'f':
            res.append('seta '+s.name+' '+parent+' '+alist)                       
            res.append('comp '+s.name+' do')                       
            res.append('comp '+s.name+' do')                       
            # write component of line element faces
        return res

    def write_gmsh(s):
        # this returns a list of lines or pts that make the component in
        # gmsh
        res = []
        
        # set the type of component
        parent = 'l'
        if isinstance(s.items[0], Line) or isinstance(s.items[0], Arc):
            parent = 'l'
        elif isinstance(s.items[0], Point):
            parent = 'p'
        
        # pull the list of entities
        alist = ','.join([str(a.id) for a in s.items])
        
        line = ''
        if parent == 'l':
            line = "Physical Line('"+s.name+"') = {" + alist + '};'
        elif parent == 'p':
            line = "Physical Point('"+s.name+"') = {" + alist + '};'
        res.append(line)
        return res

class Load_linear(Idobj):
    def __init__(s, ltype, comp, const, mult, xo, axis='x'):
        # this stores a load on a component
        # press, fx, fy, fz, ux, uy, uz
        s.ltype = ltype
        s.comp = comp
        s.const = const
        s.mult = mult
        s.xo = xo
        s.axis = axis
        Idobj.__init__(s)
    def get_press(s, x):
        res = s.const + (s.xo - x)*s.mult
        return res

    def get_list(s):
        # this returns a list of [item, val]
        res = []
        faces = s.comp.get_children()
        for f in faces:
            xs = [getattr(n, s.axis) for n in f.nodes]
            xavg = sum(xs)/len(xs)
            pval = s.get_press(xavg)
            res.append([f,pval])
        return res

    def ccx(s):
        # returns a list of lines writing the load to ccx
        res = []
        res.append('** Fluid pressure on component: '+s.comp.name)
        faces = s.comp.get_children()
        for f in faces:
            res.append('*DLOAD')
            xs = [getattr(n, s.axis) for n in f.nodes]
            xavg = sum(xs)/len(xs)
            pval = s.get_press(xavg)
            myline = '%i,P%i,%f' % (f.element.id, f.id, pval)
            res.append(myline)    
        return res

class Load(Idobj):
    def __init__(s, ltype, comp, val):
        # this stores a load on a component
        # ltype = press, fx, fy, fz, ux, uy, uz, thickness, matl
        s.ltype = ltype
        s.comp = comp
        s.val = val
        Idobj.__init__(s)

    def get_list(s):
        # this returns a list of [item, val]
        res = []
        if s.ltype == 'press':
            # pressure
            faces = s.comp.get_children()
            for f in faces:
                res.append([f,s.val])
        return res
        
    def ccx(s):
        # returns a list of lines writing the load to ccx
        res = []
        if s.ltype == 'press':
            # pressure
            res.append('** Pressure on component: '+s.comp.name)
            res.append('*DLOAD')
            faces = s.comp.get_children()
            for f in faces:
                myline = '%i,P%i,%f' % (f.element.id, f.id, s.val)
                res.append(myline)
        elif s.ltype == 'gravity':
            # gravity
            cname = s.comp.name
            res.append('** Gravity on component: '+cname)
            res.append('*DLOAD')
            res.append('%s,GRAV,%f,-1.,0.,0.' % (cname, s.val))            
        elif s.ltype[0] == 'u':
            # displacement
            axis = s.ltype[1]
            anum = {'x':1,'y':2,'z':3}
            axis = anum[axis]
            res.append('*BOUNDARY')
            res.append('%s,%i,%i,%f' % (s.comp.name, axis, axis, s.val))
        elif s.ltype[0] == 'f':
            # force
            res.append('** Force on component: '+s.comp.name)
            nodes = len(s.comp.items)
            fval = s.val/nodes
            axis = s.ltype[1]
            anum = {'x':1,'y':2,'z':3}
            axis = anum[axis]
            res.append('*CLOAD')
            res.append('%s,%i,%f' % (s.comp.name, axis, fval))
        elif s.ltype == 'thickness':
            res.append('*NODAL THICKNESS')
            res.append('%s, %f' % (s.comp.name, s.val))
        elif s.ltype == 'radps':
            radsq = s.val**2 # square the speed in radians per sec
            res.append('** Rotation (radians, value=radians^2) on component: '+s.comp.name)
            res.append('*DLOAD')
            res.append('%s,CENTRIF,%f,0.,0.,0.,0.,1.,0.' % (s.comp.name, radsq))
        elif s.ltype == 'rpm':
            rad = 2.0*pi*(s.val/60.0)
            radsq = rad**2 # square the speed in radians per sec
            res.append('** Rotation (rpm to radians, value=radians^2) on component: '+s.comp.name)
            res.append('*DLOAD')
            res.append('%s,CENTRIF,%f,0.,0.,0.,0.,1.,0.' % (s.comp.name, radsq))
        elif s.ltype == 'matl':
            mat = s.val.name
            res.append('*SOLID SECTION,ELSET='+s.comp.name+',MATERIAL='+mat)
            
        return res

    def write_cgx(s):
        # this returns a list of lines that apply loads to components in the
        # preprocessor
        res = []
        char = s.ltype[0]
        if char == 'p':
            # apply a pressure
            res.append('send '+s.comp.name+' abq pres '+str(s.val))
            fname = s.comp.name+'.dlo'
            s.set_file(fname)
        elif char == 'f':
            # aply a force
            axis = s.ltype[1]
            anum = {'x':1,'y':2,'z':3}
            ind = anum[axis]
            forces = ['0.0', '0.0', '0.0']
            forces[ind] = str(s.val)
            forces = ' '.join(forces)
            axis = anum[axis]
            res.append('send '+s.comp.name+' abq force '+forces)
            fname = s.comp.name+'.frc'
            s.set_file(fname)
        elif char == 'u':
            # apply a displacement
            axis = s.ltype[1]
            anum = {'x':'1','y':'2','z':'3'}
            axis = anum[axis]
            res.append('send '+s.comp.name+' abq spc '+axis+' '+str(s.val))
            fname = s.comp.name+'_'+axis+str(s.val)+'.bou'
            s.set_file(fname)
        return res
    

def get_ext_list(elements):
    # this returns a list of tuples of (x,y) going ccw around the selected elements
    f_bynode = {}
    begind = 1
    for e in elements:
        for f in e.faces():
            if f.ext == True:
                f_bynode[f.nodes[begind].id] = f
    
    # store results
    res = []

    # get first node
    n1 = list(f_bynode.values())
    n1 = n1[0].nodes[begind]
    res.append( (n1.y, n1.x) )
    n = n1
    
    # while the next node is not equal to the first node keep storing results
    while True:
        n = f_bynode[n.id].nodes[int(not begind)]
        if n == n1:
            break
        res.append( (n.y, n.x) )
    return res

# this class is used to make and store parts (closed areas)
class Part(Idobj):
    def __init__(s, parent):
        s.p = parent
        s.cursor = Point(0,0)
        s.matl = None
        thickness = None
        s.areas = [] # top area is the buffer
        # make the buffer
        area = s.p.areas.append( Area(s, []) )
        s.areas.append(area)
        Idobj.__init__(s)             
        s.left = None
        s.right = None
        s.top = None
        s.bottom = None
        s.file = None
        s.meshmult = 1
        s.loads = {}
    def get_item(s, item):
        # return item from part given string identifying item
        if item in ['left','right','top','bottom']:
            items = getattr(s, item)
            return items
        elif item[0] == 'P':
            # get point
            items = s.get_points()
            num = int(item[1:])
            res = [a for a in items if a.id == num]
            return res[0]
        elif item[0] == 'L':
            # get line
            items = s.get_lines()
            num = int(item[1:])
            items = [a for a in items if a.id == num]
            return items[0]
        elif item[0] == 'A':
            # get area
            items = s.areas
            num = int(item[1:])
            items = [a for a in items if a.id == num]
            return items[0]
        else:
            print('Unknown item! Please pass the name of a point, line or area!')
    def get_name(s):
        # returns part name
        return 'PART'+str(s.id)

    def set_side(s, loc, ind, axis):
        # this sets the part.left to a list of the left lines
        # loc = 'left', ind = 0, axis = 'y'
        points = s.get_points()
        # sort the points low to high
        points = sorted(points, key=lambda pt: getattr(pt, axis))
        # store the value
        val = getattr(points[ind], axis)
        res = []
        lines = s.get_lines()
        for line in lines:
            if isinstance(line, Line):
                if getattr(line.pt(0),axis) == val and getattr(line.pt(1),axis) == val:
                    # line is on the left side
                    res.append(line)
        setattr(s, loc, res)
            
    def goto(s, x, y):
        # adds a point to fea model points if it doesn't exist in part
        # otherwise just set cursor to it and return the cursor
        [pnew, already_exists] = s.make_get_pt(x, y)

        if already_exists:
            if s.areas[-1].closed == True:
                # make a new area if the old area is already closed and we're
                # goint to an existing point
                a = s.p.areas.append( Area(s, []) )
                s.areas.append( a )
        
        # return cursor
        s.cursor = pnew
        return s.cursor
    
    def draw_arc(s, x2, y2, xc, yc):
        # pass in the arc end point, arc center
        # current point is point1 in arc
        pold = s.cursor
        # make arc center point
        ctr = Point(xc,yc)
        if ctr in s.p.points:
            ctr = next(x for x in s.p.points if x == ctr)
        else:
            s.p.points.append( ctr )
        # make arc end point
        s.cursor = s.p.points.append( Point(x2,y2) )
        # make arc
        arc = s.p.lines.append( Arc(pold, s.cursor, ctr) )
        s.areas[-1].add(arc)
        return [arc, pold, s.cursor]

    def draw_line_delta(s, dx, dy):
        # draws a line given delta movements
        x = s.cursor.x + dx
        y = s.cursor.y + dy
        return s.draw_line_to(x, y)
    
    def draw_line_rad(s, dx_rad):
        # this draws a radial line of length dx
        return s.draw_line_delta(dx_rad, 0.0)
    
    def draw_line_ax(s, dy_ax):
        # this draws an axial line of length dy
        return s.draw_line_delta(0.0, dy_ax)
    
    def make_get_pt(s, x, y):
        # returns a point if it exists or a new point if it doesn't
        # return (point, alread_existed)

        pnew = Point(x,y)
        pts = s.get_points()
        pexists = False
        for p in pts:
            dist = p - pnew
            dist = dist.length()
            if dist < _acc:
                # point already exists in part, use it
                pexists = True
                pnew = p
                break
        else:
            # fired when we haven't broken out of the loop, the point is new
            pnew = s.p.points.append( pnew )
        return [pnew, pexists]
    
    def make_get_line(s, lnew):
        # returns a line if it exists or a new line if it doesn't
        # return (line, alread_existed)
        
        lnew_pset = set([p.id for p in lnew._points_])
        lines = s.get_lines()
        lexists = False
        for l in lines:
            l_pset = set([p.id for p in l._points_])
            if type(lnew) == type(l) and lnew_pset == l_pset:
                # the lines are the same type and have the same points in them
                lexists = True
                if l.pt(0) == lnew.pt(0):
                    # they are the same orientation
                    lnew = l
                else:
                    # l needs to be flipped
                    lnew = l.reverse()
                break
        else:
            # fired when we haven't broken out of the loop, the line is new
            lnew = s.p.lines.append( lnew )
        
        return [lnew, lexists]
            

    def draw_line_to(s, x, y):
        # draws a line from the current cursor poistion to the new x,y location
        
        pold = s.cursor
        [s.cursor, p_already_exists] = s.make_get_pt(x,y)
        [line, l_already_exists] = s.make_get_line( Line(pold, s.cursor) )

        s.areas[-1].add(line)
        
        # check for closure of the area
        if s.areas[-1].closed:
            s.set_side('left',0,'y')
            s.set_side('right',-1,'y')
            s.set_side('top',-1,'x')
            s.set_side('bottom',0,'x')

        return [line, pold, s.cursor]

    def get_lines(s):
        # returns lines in part
        lines = []
        for area in s.areas:
            for line in area.lines:
                if line not in lines:
                    lines.append(line)
        return lines
    
    def get_points(s):
        # returns points in part
        points = []
        for area in s.areas:
            for line in area.lines:
                for pt in line._points_:
                    if pt not in points:
                        points.append(pt)
        return points
    def get_maxlength(s):
        # returns the max length between points in an part
        pts = s.get_points()
        maxlen = 0.0
        # loop through points checking dist to next point
        for ind, p1 in enumerate(pts[:-1]):
            for p2 in pts[ind:]:
                vect = p1 - p2
                dist = vect.length()
                if dist > maxlen:
                    maxlen = dist
        return maxlen

    def area_from_pt(s, pt):
        # returns the area that the point is inside
        for area in s.areas:
            if area.pt_inside(pt):
                return area
        return None

    def fillet_lines(s, l1, l2, radius):
        # this function fillets lines in a part
        # check if the lines are touching
        tmp = s.cursor
        
        if l1.touches(l2):
            # offset the lines, assuming area is being traced clockwise
            # get the intersection point
            magnitude = radius
            l1_off = l1.offset(magnitude)
            l2_off = l2.offset(magnitude)
            ctrpt = l1_off.intersects(l2_off)
            if type(ctrpt) == type(None):
                # flip the offset deirection if lines don't intersect
                magnitude = -radius
                l1_off = l1.offset(magnitude)
                l2_off = l2.offset(magnitude)
                ctrpt = l1_off.intersects(l2_off)  
                
            # now we have an intersecting point
            print('Arc center pt is: ',ctrpt)
            p1_new = l1.arc_tang_intersection(ctrpt, magnitude)
            p2_new = l2.arc_tang_intersection(ctrpt, magnitude)
            rempt = l1.pt(1)
            
            # delete rempoint and make + store new points for the arc
            s.p.points.remove(rempt)
            [p1_new,b] = s.make_get_pt( p1_new.x, p1_new.y )
            [ctrpt,b] = s.make_get_pt( ctrpt.x, ctrpt.y )
            [p2_new,b] = s.make_get_pt( p2_new.x, p2_new.y )
            
            # make the new arc
            arc = s.p.lines.append( Arc(p1_new, p2_new, ctrpt) )
                        
            # edit the adjacent lines to replace the removed pt
            l1.set_pt(1,arc.pt(0))
            l2.set_pt(0,arc.pt(1))
            
            # put the arc in the right location in the area
            for area in s.areas:
                if l1 in area.lines:
                    area.line_insert(l1, arc)
                    print ('Arc inserted into area %i' % (area.id,))
            
            # remove point we're not using anymore
            s.p.points.remove(rempt)
            
        else:
            print ('Cannot fillet! Lines must touch!')
        # reset the cursor to where it should be
        s.cursor = tmp
        
    def point_parents(s, pt):
        # returns lines that contain a point (end point)
        res = []
        pts = s.get_points()
        lines = s.get_lines()
        for line in lines:
            if pt in line._points_:
                res.append(line)
        return res
    
    def point_line(s, pt):
        # returns the line that the point is on
        res = None
        lines = s.get_lines()
        for line in lines:
            if line.coincident(pt) == True:
                return line
        return res
    
    def line_parents(s, line):
        # returns the area parents of a line
        res = []
        for area in s.areas:
            if line in area.lines:
                res.append(area)
        return res
    
    def cut_line(s, pt, line):
        # splits a line, and updates the areas that reference the line
        # make new point real
        p = Point(pt.x,pt.y)
        pnew = None
        if p in s.p.points:
            pnew = [pt for pt in s.p.points if pt == p][0]
        else:
            pnew = s.p.points.append( p )
        # below is code to split line, requires defined point and 
        pend = line.pt(1)
        line.set_pt(1, pnew)
        lnew = s.p.lines.append( Line(pnew, pend) ) 
        areas = s.line_parents(line)
        # insert the new line into existing areas
        for area in areas:
            area.line_insert(line, lnew)
        return [pnew, lnew]

    def cut(s, pt, cvect):
        # this cuts the part, through one or more areas
        cvect.make_unit()
        
        # is the point an end point of a line?
        lines_excluded = s.point_parents(pt)
        
        # if the pt is on an existing line, add that line to the excl list
        if len(lines_excluded) == 0:
            L1 = s.point_line(pt)
            if L1 != None:
                lines_excluded.append(L1)

        # this stores the list of points to make cut new lines
        linepts = [pt]
        
        if len(lines_excluded) == 0:
            print('The given point is not on any lines, all lines will be checked for intersections!')
        elif len(lines_excluded) == 1:
            print('Cut or the origin line required')
            # add the new point and split the origin line
            line_origin = lines_excluded[0]
            [pnew, lnew] = s.cut_line(pt, line_origin)
            # set the first point of our cut line to be this actual point
            linepts[0] = pnew
        
        # make the tool cut line
        vsize = s.get_maxlength()
        endpt = pt + cvect*vsize
        cutline = Line(pt,endpt)
        
        # find the intersections
        overlaps = []
        lines = s.get_lines()
        for l in lines:
            if l not in lines_excluded:
                if isinstance(l, Line):
                    newpt = l.intersects(cutline)
                    if isinstance(newpt, Point):
                        areas = len(s.line_parents(l))
                        
                        # check if new point already exists
                        if (newpt - l.pt(0)).length() < _acc:
                            newpt = l.pt(0)
                            l = None
                        elif (newpt - l.pt(1)).length() < _acc:
                            newpt = l.pt(1)
                            l = None

                        dist = newpt - pt
                        dist = dist.length()

                        # store itersections
                        # pt, line, distance, numareas
                        d = {'pt':newpt,'line':l,'dist':dist,'areas':areas}
                        overlaps.append(d)
        
        # sort cuts by distance
        overlaps = sorted(overlaps, key=lambda k: k['dist'])
        
        # loop through them drawing the cut line
        # and closing the old and new areas
        for cut in overlaps:
            # splits the line to be cut, updates all parent areas
            pnew = None
            if cut['line'] == None:
                # don't need to cut a line at the end
                pnew = cut['pt']
            else:
                [pnew, lnew] = s.cut_line(cut['pt'], cut['line'])
            # adds the intersection point to the list of cut line points
            linepts.append(pnew)
            p1 = linepts[-2]
            p2 = linepts[-1]
            pav = p1+p2
            pav = pav*0.5
            
            # get area to cut
            area = s.area_from_pt(pav)
            # make the cut line
            cut_line = s.p.lines.append( Line(p1, p2) )
            cut1 = cut2 = cut_line
            lpre1 = area.line_from_startpt(p1)
            lpre2 = area.line_from_startpt(p2)
            
            # select the line portions that define the areas
            lind1 = area.lines.index(lpre1) 
            lind2 = area.lines.index(lpre2)
            low = min(lind1, lind2)
            high = max(lind1, lind2)
            
            # lists of lines for areas
            beg = area.lines[:low]
            end = area.lines[high:]
            mid = area.lines[low:high]
            
            # reverse lines if needed
            p1 = beg[-1].pt(1)
            if cut1.pt(0) != p1:
                cut1 = cut1.reverse()
            p1 = mid[-1].pt(1)
            if cut2.pt(0) != p1:
                cut2 = cut2.reverse()
            
            alist_curr = beg + [cut1] + end
            alist_other = mid + [cut2]
            
            area.update(alist_curr)
            anew = s.p.areas.append( Area(s, alist_other) )
            s.areas.append(anew)

            if cut['areas'] == 1:
                # break out of for loop, we hit a part edge
                break
        
    def chunk_area(s, area):
        # this chunks an area, cutting it into smaller meshable areas
        # walk around the area looking at each point        

        # store the cuts first, then cut after
        cuts = [] # each item is a dict with a pt and vect in it
        for ind, line in enumerate(area.lines):
            L1 = area.lines[ind-1]
            L2 = line
            pt = L1.pt(1)
            perp1 = L1.get_abc_perp(pt)
            perp2 = L2.get_abc_perp(pt)
            v1 = L1.get_perp_vec(pt)
            v2 = L2.get_perp_vec(pt)
            # flip these vectors later to make them cut the area(s)
            ang = v1.ang_bet_deg(v2)
            case = None
            d = {}
            if ang == 0.0:
                # cut this
                case = 'tangent'
                d = {'pt':pt, 'vect':v1*-1}
                cuts.append(d)
            elif ang > 0:
                # cut this
                case = 'internal_corner'
                d = {'pt':pt, 'vect':v1*-1}
                cuts.append(d)
                d = {'pt':pt, 'vect':v2*-1}
                cuts.append(d)
            elif ang < 0:
                case = 'external_corner'
                # do not split these
        
        # do the cuts
        for cut in cuts:
            print('--------------------')
            print('Cut pt:',cut['pt'])
            print('Cut vect: ',cut['vect'])
            s.cut(cut['pt'], cut['vect'])
            
    def chunk(s):
        # this method splits the part area into 3-5 sided area pieces
        for area in s.areas:
            if area.closed:
                if len(area.lines) > 5:
                    # we need to chunk this area, it has too many lines to mesh
                    s.chunk_area(area)
                else:
                    print('Area %i was not chunked because it had < 6 sides.'%(area.id))
        # store the left, right, top, and bottom lines        
        s.set_side('left',0,'y')
        s.set_side('right',-1,'y')
        s.set_side('top',-1,'x')
        s.set_side('bottom',0,'x')
    def set_matl(s, mat):
        # this sets the matl of the part
        s.matl = mat

    def plot_pressures(s, fname='', display=True):
        # plot the pressures on elements

        pts = s.get_points()
        radials = [p.x for p in pts]
        axials = [p.y for p in pts]
        
        # plot all elements and store length, length determins min pressure arrow
        if hasattr(s, 'elements'):
            
            # plotting elements
            fig = plt.figure()
            ax = fig.add_subplot(111,aspect='equal')

            patches = []
            face_len = []
            for e in s.elements:
                face_len.append( e.face[1].length() )
                xycoords = e.get_corner_nodes()
                polygon = Polygon(xycoords, closed=True)
                patches.append(polygon)
            p = PatchCollection(patches, edgecolors=_ecolor, facecolors=_fcolor)
            ax.add_collection(p)
                        
            # average face length is min arrow length
            face_len = sum(face_len)/len(face_len)            

            # store pressures we'll want to plot
            # this is a list of lists [face, pval]
            plist = []
            for load in s.p.loads[s.p.time]:
                if load.ltype in  ['press','press_fluid']:
                    plist += load.get_list()
            pressures = [pval for [face,pval] in plist]

            # check min and max bounds
            pmin = min(pressures)
            pmax = max(pressures)
            mult = face_len/abs(pmin)    # mult to go from pressure to length
            
            # make tick list for later plot, and color map
            tick_list = [pmin]  # default to plot one val
            cmap = colors.ListedColormap(['b', 'b']) # default to plot one val
            if pmax != pmin:
                # we have a range of values we're plotting
                tick_list = frange(pmin,pmax,(pmax-pmin)/8)
                cmap = plt.cm.jet
                        
            # set color contours for arrows
            cNorm  = colors.Normalize(vmin=pmin, vmax=pmax)
            scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=cmap)            
            scalarMap._A = [] # need to set this for it to work
            
            # make arrows
            for [face, pval] in plist:
                [p1, unit] = face.get_mnorm()
                pdelta = unit*(mult*abs(pval))
                p2 = p1 + pdelta
                
                # assuming positive pressure, arrow points to face p2->p1
                pstart = p2
                delta = p1 - p2
                if pval < 0:
                    pstart = p1
                    delta = p2 - p1                    
                radials.append(p2.x)
                axials.append(p2.y)

                hw = face_len*0.2
                hl = face_len*0.3
                colorVal = scalarMap.to_rgba(pval)
                plt.arrow(pstart.y,  #x1
                          pstart.x,  # y1
                          delta.y, # x2 - x1
                          delta.x, # y2 - y1
                          color=colorVal, head_width=hw, head_length=hl,
                          length_includes_head=True)                
                            
            # set the horizontal and vertical axes
            vert = max(radials) - min(radials)
            horiz = max(axials) - min(axials)
            vadder = (vert)/5
            hadder = (horiz)/5    
            (vmax,vmin) = (max(radials)+vadder, min(radials)-vadder)
            (hmax,hmin) = (max(axials)+hadder, min(axials)-hadder)

            # set units
            [d_unit, p_unit, t_unit] = s.p.get_units('dist', 'stress', 'time')

            # set plot axes
            tstr = '%s pressures %s\nTime=%f%s' % (s.get_name(), p_unit, s.p.time,
                                                 t_unit)
            plt.title(tstr)
            plt.xlabel('axial, y'+d_unit)
            plt.ylabel('radial, x'+d_unit)

            # set min and max vertical and axial limits
            plt.xlim(hmin, hmax)
            plt.ylim(vmin, vmax)            
            
            # set the colorbar
            cbar = plt.colorbar(scalarMap, orientation='vertical', ticks=tick_list)
            
            if fname != '':
                # save the image
                fname += '.png'
                if _dpi != None:
                    plt.savefig(fname, dpi=_dpi, bbox_inches='tight')
                else:
                    plt.savefig(fname, bbox_inches='tight')
            
            if display:
                plt.tight_layout()
                plt.show()        

            # remove all figures
            plt.close()

        else:
            # part has not been meshed yet
            res = 'Part: %s does not have any elements! ' % (s)
            res += 'Try meshing it with model.mesh(1)' 
            print(res)

        
    def plot_geometry(s, fname='', display=True):
        # this method plots the part
        # check out: http://nickcharlton.net/posts/drawing-animating-shapes-matplotlib.html
        
        pts = s.get_points()
        nax=[pt.y for pt in pts]
        nrad=[pt.x for pt in pts]
        n_id=[pt.get_name() for pt in pts]
        
        # need to go through lines making list of points, don't have to be
        # sequential

        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        #fig, ax = plt.subplots()
        ax.scatter(nax, nrad)
        
        # plot points
        for i, txt in enumerate(n_id):
            ax.annotate(txt, (nax[i],nrad[i]))
        
        #plot lines
        lines = s.get_lines()
        for l in lines:
            l.plot(ax)
        
        #plot area ids
        for a in s.areas:
            if a.closed == True:
                pax = a.centroid.y
                prad = a.centroid.x
                ax.annotate(a.get_name(), (pax,prad))

        # set the horizontal and vertical axes
        vert = max(nrad) - min(nrad)
        horiz = max(nax) - min(nax)
        vadder = (vert)/5
        hadder = (horiz)/5
        (vmax,vmin) = (max(nrad)+vadder, min(nrad)-vadder)
        (hmax,hmin) = (max(nax)+hadder, min(nax)-hadder)
        plt.xlim(hmin, hmax)
        plt.ylim(vmin, vmax)

        # set units
        [d_unit] = s.p.get_units('dist')            
        
        # show plot
        plt.title(s.get_name()+' geometry')
        plt.xlabel('axial, y'+d_unit)
        plt.ylabel('radial, x'+d_unit)
        ax.set_aspect('equal')

        if fname != '':
            # save the image
            fname += '.png'
            if _dpi != None:
                plt.savefig(fname, dpi=_dpi, bbox_inches='tight')
            else:
                plt.savefig(fname, bbox_inches='tight')
        
        if display:
            plt.tight_layout()
            plt.show()
        
        # remove all figures
        plt.close()

    def plot_elements(s, fname='', display=True):
        # plots the elements in the part

        pts = s.get_points()
        radials = [p.x for p in pts]
        axials = [p.y for p in pts]
        vadder = (max(radials) - min(radials))/5
        hadder = (max(axials) - min(axials))/5

        (vmax,vmin) = (max(radials)+vadder, min(radials)-vadder)
        (hmax,hmin) = (max(axials)+hadder, min(axials)-hadder)
        
        if hasattr(s, 'elements'):
            # plotting elements
            fig, ax = plt.subplots()            
            patches = []
            for e in s.elements:
                xycoords = e.get_corner_nodes()
                polygon = Polygon(xycoords, closed=True)
                patches.append(polygon)
            p = PatchCollection(patches, facecolors=_fcolor, edgecolors=_ecolor)
            ax.add_collection(p)
            
            # set units
            [d_unit] = s.p.get_units('dist')            
            
            plt.title(s.get_name()+' elements')
            plt.xlabel('axial, y'+d_unit)
            plt.ylabel('radial, x'+d_unit)
            plt.axis('scaled')
            plt.xlim(hmin, hmax)
            plt.ylim(vmin, vmax)
            
            if fname != '':
                # save the image
                fname += '.png'
                if _dpi != None:
                    plt.savefig(fname, dpi=_dpi, bbox_inches='tight')
                else:
                    plt.savefig(fname, bbox_inches='tight')
            
            if display:
                plt.tight_layout()
                plt.show()        

            # remove all figures
            plt.close()

        else:
            # part has not been meshed yet
            res = 'Part: %s does not have any elements! ' % (s)
            res += 'Try meshing it with model.mesh(1)' 
            print(res)
            
    def __str__(s):
        # prints a string representation of the part
        val = 'Part, id=%i name=%s' % (s.id, s.get_name())
        return val

# stores a node
class Node(object):
    def __init__(s, node_num, x, y, z):
        s.id = node_num
        s.x = x
        s.y = y
        s.z = z
        s.order = 1
    def set_order(s, val):
        # 1 or 2
        s.order = val
    def __hash__(s):
        return s.id
    def __eq__(s, other):
        if isinstance(other, Node):
            return s.id == other.id
        else:
            return False
    def ccx(s):
        # writes out line defining node for ccx inp writing
        val = '%i, %f, %f, %f' % (s.id, s.x, s.y, s.z)
        return val
    def __str__(s):
        val = 'Node, id %i, (x,y)=(%f,%f)' % (s.id, s.x,s.y)
        return val


#stores an element face
class Face(object):
    def __init__(s, fnum, n1, n2, element):
        s.id = fnum
        s.element = element
        s.nodes = [n1,n2]
        s.nmid = None
        s.ext = False   # stores if face is external
    def length(s):
        # return the length of the face
        p0 = Point(s.nodes[0].x,  s.nodes[0].y)
        p1 = Point(s.nodes[1].x,  s.nodes[1].y)
        vect = p1 - p0
        return vect.length()
    def get_mnorm(s):
        # return midpt, normal
        p0 = Point(s.nodes[0].x,  s.nodes[0].y)
        p1 = Point(s.nodes[1].x,  s.nodes[1].y)
        midpt = p0 + p1
        midpt = midpt*0.5
        norm = p1 - p0
        norm.rot_ccw_deg(90)
        norm.make_unit()
        return [midpt, norm]
    def set_ext(s):
        # set the face to be external
        s.ext = True
    def set_nmid(s, nmid):
        s.nmid = nmid
    def __eq__(s, other):
        # check equality to other faces
        my_nodes = set([s.nodes[0].id,s.nodes[1].id])
        other_nodes = set([other.nodes[0].id,other.nodes[1].id])
        if my_nodes == other_nodes:
            return True
        else:
            return False
    def __str__(s):
        n1 = s.nodes[0].id
        n2 = s.nodes[1].id
        fid = s.id
        eid = s.element.id
        s = 'Face %i on element %i with nodes [%i,%i]' % (fid, eid, n1, n2)
        return s

# stores an element
class Element(object):
    def __init__(s, enum, etype, nlist):
        s.id = enum
        s.etype = etype
        s.node = {}
        s.face = {}

        # figure out the number of faces
        fnum = len(nlist)
        if fnum > 4:
            fnum = int(fnum/2)

        # store nodes
        for (ind, node) in enumerate(nlist):
            s.node[ind+1] = node
            if ind >= fnum:
                node.set_order(2)
        
        # store faces
        n = nlist[0:fnum]+[nlist[0]]
        for i in range(fnum):
            n1 = n[i]
            n2 = n[i+1]
            f = Face(i+1, n1, n2, s)
            s.face[i+1] = f
    def get_tris(s):
        # this returns the triangles for plotting of results
        # triangles must be closed in a CCW direction
        res = []
        numnodes = len(s.nodes()) 
        if numnodes == 3:
            # triangle, first order = 1 triangle
            res.append( [s.node[1].id, s.node[2].id, s.node[3].id] )
        elif numnodes == 4:
            # quad, first order = 2 trianges
            res.append( [s.node[1].id, s.node[2].id, s.node[3].id] )
            res.append( [s.node[1].id, s.node[3].id, s.node[4].id] )
        elif numnodes == 6:
            # tri, second order = 4 triangles
            res.append( [s.node[1].id, s.node[4].id, s.node[6].id] )
            res.append( [s.node[4].id, s.node[2].id, s.node[5].id] )
            res.append( [s.node[6].id, s.node[5].id, s.node[3].id] )
            res.append( [s.node[6].id, s.node[4].id, s.node[5].id] )
        elif numnodes == 8:
            # quad, second order = 6 triangles
            res.append( [s.node[1].id, s.node[5].id, s.node[8].id] )
            res.append( [s.node[5].id, s.node[2].id, s.node[6].id] )
            res.append( [s.node[6].id, s.node[3].id, s.node[7].id] )
            res.append( [s.node[8].id, s.node[7].id, s.node[4].id] )
            res.append( [s.node[8].id, s.node[5].id, s.node[6].id] )
            res.append( [s.node[8].id, s.node[6].id, s.node[7].id] )
        return res
    
    def get_area(s):
        # returns the area of the element
        asum = 0.0
        for face in s.faces():
            n1 = face.nodes[0]
            n2 = face.nodes[1]
            asum += n2.y*n1.x - n2.x*n1.y
        asum = asum*0.5
        return asum
    
    def faces(s):
        #returns element faces
        return s.face.values()
    
    def nodes(s):
        # return element nodes
        return s.node.values()
    
    def get_corner_nodes(s):
        # this is used for plotting so reverse x and y
        res = []
        nnum = len(s.faces())
        for ind in range(1,nnum+1):
            node = s.node[ind]
            res.append( [node.y, node.x] )
        return res
    
    def set_etype(s, etype):
        # changes the element type to passes axisym, plstress, or plstrain
        fnum = len(s.faces())
        shape = 'tri'
        if fnum == 4:
            shape = 'quad'
        order = '2'
        if len(s.nodes()) == fnum:
            order = '1'
        estr = shape+order+etype
        s.etype = _ccx_elements[estr]
        
    def ccx(s):
        # returns ccx string defining element
        nids = [str(n.id) for n in s.node.values()]
        val = str(s.id)+', '+', '.join(nids)
        return val

class Results_File(object):
    def __init__(s, solved_model, fname):
        # this stores a results file in python
        s.p = solved_model
        s.fname = fname
        s.steps = [] # this stores a list of time steps
        s.results = {} # stores results, nested dicts
        s.time = -1
        s.read_frd()
        #s.read_dat() # triggers read
    def nplot(s, field, fname='', display=True, levels=21, gradient=False, gmult=1.0):
        # plot the results in the given field for the selected object

        # make a list of the selected nodes and elements
        sel = {'nodes':set( [] ),'elements':[]}
        selected = s.p.p.selected

        # store selected nodes and elements
        for item in selected:
            if isinstance(item, Element):
                sel['elements'] += [item]
                sel['nodes'].update( item.nodes() )
            elif isinstance(item, Part) or isinstance(item, Area):
                sel['elements'] += item.elements
                sel['nodes'].update(item.nodes)
        
        # sort nodes low to high so index is correct
        # WITH SUBSET I'LL NEED TO REWRITE THIS WITH A ID TO INDEX MAPPING
        sel['nodes'] = list(sel['nodes'])
        sel['nodes'] = sorted(sel['nodes'], key=lambda k: k.id)
        
        # store results at nodes
        axials = []
        radials = []
        zs = []
        id_to_ind = {}
        for n in sel['nodes']:
            id_to_ind[n.id] = len(axials)
            ax = n.y + gmult*s.results[s.time]['node'][n.id]['uy']
            rad = n.x + gmult*s.results[s.time]['node'][n.id]['ux']
            axials.append(ax)
            radials.append(rad)
            zs.append(s.results[s.time]['node'][n.id][field])

        # make a list of triangles, given by indices, looping anticlockwise
        triangles = []
        for e in sel['elements']:
            tris = e.get_tris()     # list of triangle nodes defined by node id
            for t in tris:
                for ind, nid in enumerate(t):
                    t[ind] = id_to_ind[nid]     # convert id to index
            triangles += tris
        
        # check to see if selected nodes and elements are
        # in the parent model's nodes and elements
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # need to set tick list here
        vmin = min(zs)
        vmax = max(zs)
        tick_list = [vmin]
        if vmax != vmin:
            # we have a range of values we're plotting
            tick_list = frange(vmin,vmax,(vmax-vmin)/levels)        
        
        # plot using a gradient(shaded) or levels
        if gradient:
            # This one is shaded
            plt.tripcolor(axials, radials, triangles, zs, shading='gouraud')
        else:
            # this one is not shaded
            plt.tricontourf(axials, radials, triangles, zs, levels=tick_list)

        # code required for the colorbar
        cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = colors.ListedColormap(['b', 'b']) # default to plot one val
        if vmax != vmin:
            # we have a range of values we're plotting
            if gradient:
                cmap = plt.cm.jet
            else:
                cmap = plt.get_cmap('jet', levels)
        scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=cmap)            
        scalarMap._A = [] # need to set this for it to work
        plt.colorbar(scalarMap, orientation='vertical', ticks=tick_list)
            
        # set the horizontal and vertical axes
        vert = max(radials) - min(radials)
        horiz = max(axials) - min(axials)
        vadder = (vert)/5
        hadder = (horiz)/5    
        (vmax,vmin) = (max(radials)+vadder, min(radials)-vadder)
        (hmax,hmin) = (max(axials)+hadder, min(axials)-hadder)
        plt.xlim(hmin, hmax)
        plt.ylim(vmin, vmax)
        
        # set units
        [f_unit, d_unit, t_unit] = s.p.p.get_units(field, 'dist', 'time')

        # set plot axes
        plt.title('%s%s\nTime=%f%s' % (field, f_unit, s.time, t_unit))
        plt.xlabel('axial, y'+d_unit)
        plt.ylabel('radial, x'+d_unit)
        ax.set_aspect('equal')
        if gmult != 1:
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])

        if fname != '':
            # save the image
            fname += '.png'
            if _dpi != None:
                plt.savefig(fname, dpi=_dpi, bbox_inches='tight')
            else:
                plt.savefig(fname, bbox_inches='tight')
            print('File %s was saved' % (fname))
        
        if display:
            plt.tight_layout()
            plt.show()
            
        # remove all figures
        plt.close()
        
    def set_time(s, val):
        # sets the current time to time val
        s.time = val
        print('Results file time set to: %f' % (s.time))

    def utot(s, vals):
        # computes sum of the squares
        res = [a**2 for a in vals]
        res = ( sum(res) )**0.5
        return res
    
    def seqv(s, vals):
        # calculates the equivalent stress
        [s11,s22,s33,s12,s13,s23] = vals
        a = s11 - s22
        b = s22 - s33
        c = s33 - s11
        d = s12**2 + s23**2 +s13**2
        res = (0.5*(a**2 + b**2 + c**2 +6*d))**0.5
        return res

    def principals(s, vals):
        # calculates and returns principal stresses, S1, S2, S3
        [s11,s22,s33,s12,s13,s23] = vals
        a = 1
        b = (s11 + s22 + s33)*-1.0
        c = (s11*s22 + s11*s33 + s22*s33 - s12**2 - s13**2 - s23**2)
        d = (s11*s22*s33 + 2*s12*s13*s23 - s11*(s23**2) - s22*(s13**2) - s33*(s12**2))*-1.0
        res = list(roots([a,b,c,d]))
        res = sorted(res, reverse=True)
        
        '''
        [I1, I2, I3] = [b,c,d]
        phi = (1.0/3)*acos( (2*(I1**3)-9*I1*I2+27*I3) / (2*(I1**2 - 3*I2)**(1.5)) )
        s1 = I1/3 + (2/3)*sqrt(I1**2 - 3*I2)*cos(phi)
        s2 = I1/3 + (2/3)*sqrt(I1**2 - 3*I2)*cos(phi + 2*pi/3)
        s3 = I1/3 + (2/3)*sqrt(I1**2 - 3*I2)*cos(phi + 4*pi/3)
        res2 = [-s1, -s2, -s3]
        res2 = sorted(res2, reverse=True)
        '''
        return res
    
    def get_nmax(s, field):
        # returns the max value of a given results field
        res = [ ndict[field] for ndict in s.results[s.time]['node'].values() ]
        return max(res)

    def get_nmin(s, field):
        # returns the min value of a given results field
        res = [ ndict[field] for ndict in s.results[s.time]['node'].values() ]
        return min(res)
    
    def get_fsum(s, item):
        # returns fsum on nodes under given item (point or line)
        (fx,fy,fz) = ( [],[],[] )
        nodes = None
        if isinstance(item, Part) or isinstance(item, Part):
            nodes = item.nodes
        elif isinstance(item, Point) or isinstance(item, Line) or isinstance(item, Arc):
            nodes = item.nall

        nodes = [n.id for n in nodes]
        for n in nodes:
            x = s.results[s.time]['node'][n]['fx']
            y = s.results[s.time]['node'][n]['fy']
            z = s.results[s.time]['node'][n]['fz']
            if x != 0 or y != 0 or z != 0:
                #print('Node %i, (fx, fy, fz) = (%12.11f,%12.11f,%12.11f)' % (n, x, y, z))
                fx.append(x)
                fy.append(y)
                fz.append(z)
        fx = sum(fx)
        fy = sum(fy)
        fz = sum(fz)
        return [fx, fy, fz]
    
    def get_emax(s, stype):
        # returns the max stress of given type
        res = 0.0
        for step in s.steps:
            for (enum, edict) in s.results[step]['element'].items():
                for (ipnum, ipdict) in edict.items():
                    if ipdict[stype] > res:
                        res = ipdict[stype]
        return res

    def get_vals(s, fstr, line):
        # this returns a list of items based on an input format string
        res = []
        fstr = fstr.split(',')
        for item in fstr:
            if item[0] == "'":
                # strip off the char quaotes
                item = item[1:-1]
                # this is a string entry, grab the val out of the line
                ind = len(item)
                fwd = line[:ind]
                line = line[ind:]
                res.append(fwd)
            else:
                # format is: 1X, A66, 5E12.5, I12
                # 1X is number of spaces                
                (m,c) = (1, None)
                m_pat = re.compile(r'^\d+') # find multiplier
                c_pat = re.compile(r'[XIEA]') # find character
                if m_pat.findall(item) != []:
                    m = int(m_pat.findall(item)[0])
                c = c_pat.findall(item)[0]
                if c == 'X':
                    # we are dealing with spaces, just reduce the line size
                    line = line[m:]
                elif c == 'A':
                    # character string only, add it to results
                    fwd = line[:m].strip()
                    line = line[m:]
                    res.append(fwd)
                else:
                    # IE, split line into m pieces
                    w_pat = re.compile(r'[IE](\d+)') # find the num after char
                    w = int(w_pat.findall(item)[0])
                    for i in range(m):
                        # only add items if we have enough line to look at
                        if w <= len(line):
                            substr = line[:w]
                            line = line[w:]
                            substr = substr.strip() # remove space padding
                            if c == 'I':
                                substr = int(substr)
                            elif c == 'E':
                                substr = float(substr)
                            res.append(substr)
        return res
                
    def read_frd(s):
        # this reads a ccx frd results file which contains nodal results
        if os.path.isfile(s.fname):
            
            print('Reading results file: '+s.fname)
            f = open(s.fname)
            mode = 'off'
            time = 0.0
            FORMAT = None
            rfstr = ''
            node_ids = []            
            while True:
                line = f.readline()
                if not line:
                    break
                
                #-------------
                # set the modes
                #--------------
                if '2C' in line:
                    # we are in nodal definition
                    fstr = "1X,'   2','C',18X,I12,37X,I1"
                    t = s.get_vals(fstr, line)
                    [KEY,CODE,NUMNOD,FORMAT] = t

                    # set results format to short, long or binary
                    # only short and long are parsed so far
                    if FORMAT == 0:
                        # short format
                        rfstr = "1X,'-1',I5,3E12.5"
                    elif FORMAT == 1:
                        # long format
                        rfstr = "1X,'-1',I10,3E12.5"
                    elif format == 2:
                        # binary
                        pass
                    
                    line = f.readline()
                    mode = 'nodes'
                    print('Reading '+mode)

                elif '1PSTEP' in line:
                    # we are in a results block
                    
                    # next line has the time in it
                    line = f.readline()
                    fstr = "1X,' 100','C',6A1,E12.5,I12,20A1,I2,I5,10A1,I2"
                    t = s.get_vals(fstr, line)
                    [KEY,CODE,SETNAME,VALUE,NUMNOD,TEXT,ICTYPE,NUMSTP,ANALYS,FORMAT] = t

                    # set results format to short, long or binary
                    # only short and long are parsed so far
                    if FORMAT == 0:
                        # short format
                        rfstr = "1X,I2,I5,6E12.5"
                    elif FORMAT == 1:
                        # long format
                        rfstr = "1X,I2,I10,6E12.5"
                    elif format == 2:
                        # binary
                        pass

                    # set the time
                    time = VALUE
                    if time not in s.steps:
                        s.steps.append(time)
                    if time not in s.results:
                        s.results[time] = {'node':{},'element':{}}
                        for nid in node_ids:
                            # make dict for each node
                            s.results[time]['node'][nid] = {}
                        
                    # get the name to determine if stress or displ
                    line = f.readline()
                    fstr = "1X,I2,2X,8A1,2I5"
                    t = s.get_vals(fstr, line)
                    [KEY, NAME,NCOMPS,IRTYPE] = t
                    
                    if NAME == 'DISPR':
                        mode = 'displ'
                        # read 4 lines to get to data
                        for i in range(4):
                            f.readline()
                    elif NAME == 'STRESSR':
                        mode = 'stress'
                        # read 6 lines to get to data
                        for i in range(6):
                            f.readline()
                    elif NAME == 'TOSTRAIR':
                        mode = 'strain'
                        # read 6 lines to get to data
                        for i in range(6):
                            f.readline()
                    elif NAME == 'FORCR':
                        mode = 'force'
                        # read 4 lines to get to data
                        for i in range(4):
                            f.readline()
                            
                    print('Reading '+mode+' storing: '+','.join(_resfields[mode]))
                    line = f.readline()
                
                #----------------------------
                # Store results
                #----------------------------
                
                # reset the read mode if we hit an end of record
                if line[:3] == ' -3':
                    mode = 'off'
                
                if mode == 'nodes':
                    # node definition, store node numbers only
                    t = s.get_vals(rfstr, line)
                    [KEY, NODE, x, y, z] = t
                    node_ids.append(NODE)
                    
                elif mode == 'displ':  
                    # displacements
                    t = s.get_vals(rfstr, line)
                    [KEY, NODE, ux, uy, uz] = t                    
                    labs = _resfields[mode]
                    vals = [ux, uy, uz]
                    utot = s.utot(vals)
                    vals.append(utot)
                    #print(vals)
                    d = s.results[time]['node'][NODE]
                    for (k, val) in zip(labs, vals):
                        d[k] = val

                elif mode == 'stress':
                    # stresses
                    t = s.get_vals(rfstr, line)
                    [KEY, NODE, sx, sy, sz, sxy, syz, szx] = t
                    labs = _resfields[mode]
                    vals = [sx, sy, sz, sxy, syz, szx]
                    seqv = s.seqv(vals)
                    [s1, s2, s3] = s.principals(vals)
                    vals.append(seqv)
                    vals += [s1, s2, s3]
                    #print(vals)
                    d = s.results[time]['node'][NODE]
                    for (k, val) in zip(labs, vals):
                        d[k] = val

                elif mode == 'strain':
                    # strains
                    t = s.get_vals(rfstr, line)
                    [KEY, NODE, ex, ey, ez, exy, eyz, ezx] = t
                    labs = _resfields[mode]
                    vals = [ex, ey, ez, exy, eyz, ezx]
                    eeqv = s.seqv(vals)
                    [e1, e2, e3] = s.principals(vals)
                    vals.append(eeqv)
                    vals += [e1, e2, e3]
                    #print(vals)
                    d = s.results[time]['node'][NODE]
                    for (k, val) in zip(labs, vals):
                        d[k] = val
                
                elif mode == 'force':
                    # reaction forces
                    t = s.get_vals(rfstr, line)
                    [KEY, NODE, fx, fy, fz] = t
                    labs = _resfields[mode]
                    vals = [fx, fy, fz]
                    d = s.results[time]['node'][NODE]
                    for (k, val) in zip(labs, vals):
                        d[k] = val                    

            f.close()
            print('The following times have been read: ',s.steps)
            print('Done reading file: %s' % (s.fname))
            s.set_time(s.steps[0])
            
        else:
            # Show an error
            print("Error: %s file not found" % fname)        

    def read_dat(s):
        # reads the results file
        if os.path.isfile(s.fname):
            
            f = open(s.fname)
            readtype = 'off'
            time = 0.0
            while True:
                line = f.readline()
                if not line:
                    break
                line = line.strip()
                
                # check for read flags for displ or stress
                if 'displacement' in line or 'stress' in line:
                    words = line.split()
                    time = float(words[-1])
                    if time not in s.steps:
                        s.steps.append(time)
                    if time not in s.results:
                        s.results[time] = {'node':{},'element':{}}
                    readtype = 'displ'
                    if 'stress' in line:
                        readtype = 'stress'
                    f.readline()
                    line = f.readline().strip()
                
                # reset the read type if we hit a blank line
                if line == '':
                    readtype = 'off'
                
                # store displacement results
                if readtype == 'displ':                   
                    w = line.split()
                    labels = ['ux', 'uy', 'uz', 'utot']
                    vals = [float(a) for a in w[1:]]
                    utot = s.utot(vals)
                    vals.append(utot)
                    #print(vals)
                    node_num = int(w[0])
                    d = {}
                    for (key, val) in zip(labels, vals):
                        d[key] = val
                    s.results[time]['node'][node_num] = d

                elif readtype == 'stress':
                    # store stress results
                    w = line.split()
                    labs = ['Sx', 'Sy', 'Sz', 'Sxy', 'Sxz', 'Syz', 'Seqv']
                    vals = [float(a) for a in w[2:]]
                    seqv = s.seqv(vals)
                    vals.append(seqv)
                    #print(vals)
                    enum = int(w[0])
                    ipnum = int(w[1])
                    d = {}
                    for (key, val) in zip(labs, vals):
                        d[key] = val
                    if enum not in s.results[time]['element']:
                        s.results[time]['element'][enum] = {}
                    s.results[time]['element'][enum][ipnum] = d

            f.close()
            print('The following times have been read: ',s.steps)
            print('Results from file: %s have been read.' % (s.fname))
            
        else:
            # Show an error
            print("Error: %s file not found" % fname)        

class Model(Idobj):
    def __init__(s, parent, parts, mtype):
        # parts can be a single part or a list of parts
        # this class stores models and results, model type=
        #     struct, struct-therm, therm, cfd
        s.p = parent
        if not isinstance(parts, list):
            parts = [parts]
        s.parts = parts
        s.mtype = mtype
        s.rfile = None
        Idobj.__init__(s)
    def get_ntxt(s, nodes):
        # returns text defining all nodes
        res = []
        res.append('*NODE, NSET=Nall')
        for n in nodes:
            res.append(n.ccx())
        return res
    def get_etxt(s, elements):
        # returns text defining all elements
        res = []
        types = set([e.etype for e in elements])
        for t in types:
            tname = t
            if len(types) == 1:
                tname = 'Eall'
            res.append('*ELEMENT, TYPE='+t+', ELSET='+tname)            
            eset = [e for e in elements if e.etype == t]
            for e in eset:
                res.append(e.ccx())
        return res
    def get_ctxt(s, components):
        # returns text defining compnents of nodes or elements        
        res = []
        for c in components:
            res += c.ccx()
        return res
    def get_eset(s, name, elements):
        # returns text defining compnents of elements

        res = []
        items_per_line = 6
        res.append('*ELSET,ELSET='+name)
        grouped_els = chunk_list(elements, items_per_line)
        for group in grouped_els:
            item_ids = [str(e.id) for e in group]
            line = ', '.join(item_ids)
            if group != grouped_els[-1]:
                line += ','
            res.append(line)
        return res

    def solve(s):
        # exports a calculix file ccx input file
        inp = []
        
        # store what results we'll be outputting for eac type of analysis
        out_el = {}
        out_el['struct'] = 'E,S' # strain, stress
        out_node = {}
        out_node['struct'] = 'RF,U' # reaction forces, displacement
        
        if s.mtype == 'struct':
            # store nodes and elements
            N = []
            E = []
            C = []
            P = []

            # store all loads in the parts in this model
            these_loads = s.p.loads
            
            # store all nodes, elements, and part element sets
            for part in s.parts:
                E += part.elements
                N += part.nodes
                #P += s.get_eset(part.get_name(), part.elements)
            
            # store all nodal components
            for time in these_loads:
                for load in these_loads[time]:
                    if load.ltype not in ['press', 'press_fluid']:
                        C.append(load.comp)
            
            N = s.get_ntxt(N)
            E = s.get_etxt(E)
            C = s.get_ctxt(C)
            
            # nodes
            inp += N
            # elements
            inp += E
            # part components
            inp += P
            # load components
            inp += C
            
            # read in all materials
            for matl in s.p.matls:
                inp += matl.ccx()
            
            # write all steps and loads
            for time in these_loads:
                if time == 0:
                    # this is for thicknesses and materials
                    for load in these_loads[time]:
                        inp += load.ccx()
                else:
                    # only write times >= 1
                    inp.append('*STEP')
                    inp.append('*STATIC')
                    
                    for load in these_loads[time]:
                        inp += load.ccx()                

                    # make output frd file for cgx
                    inp.append('*EL FILE')
                    inp.append(out_el[s.mtype])
                    inp.append('*NODE FILE')
                    inp.append(out_node[s.mtype])
            
                    # end step
                    inp.append('*END STEP')
            
            '''
            # set the output files for writing
            inp.append('*NODE PRINT,NSET=Nall')
            inp.append('U')
            inp.append('*EL PRINT,ELSET=Eall')
            inp.append('S')    
            '''            
            
            # write CCX inp file to the local directory
            fname = s.p.fname+'.inp'
            f = open(fname,'w')
            for line in inp:
                #print (line)
                f.write(line+'\n')        
            print ('File: '+ fname + ' was written')
            f.close()
            
            # run file
            p = subprocess.Popen("%s %s" % (_ccx, s.p.fname))
            p.wait()
            print('Solving done!')
            
            # read the results file in
            s.rfile = Results_File(s, s.p.fname+'.frd')
            s.p.select(list(s.parts))

class FeaModel():
    def __init__(s, fname):
        s.fname = fname
        s.points = Item_List()
        s.lines = Item_List()
        s.areas = Item_List()
        s.parts = Item_List()
        s.matls = Item_List()
        s.components = Item_List()
        s.loads = {}    # store loads by time
        s.models = Item_List()
        s.nodes = ID_List()
        s.elements = ID_List()
        s.faces = []
        s.selected = None   # this selected set is what we'll plot
        s.time = 1.0 # 0 is model set-up, 1 is first step, etc
        s.units = {}
    def select(s, items):
        # set the selected item to the passed item or list of elements
        if isinstance(items, list):
            s.selected = items
        else:
            s.selected = [items]

    def set_units(s, dist='m', temp='K'):
        # sets the units that will be displayed when plotting
        keys = ['displ','force','stress','temp', 'density', 'time']
        
        m_N = ['m','N','Pa','K', 'kg/(m^3)', 's']
        mm_N = ['mm','N','MPa','K','tonne/(mm^3)', 's']
        in_lbf = ['in','lbf','psi','R', 'slinch/(in^3)', 's']
        ft_lbf = ['ft','lbf','psf','R', 'slug/(ft^3)', 's']
        dicts = [m_N, mm_N, in_lbf, ft_lbf]
        
        # pick correct dict list
        dist_unit = [d[0] for d in dicts]
        ind = dist_unit.index(dist)
        vals = dicts[ind]
        
        # set values
        d = dict( zip(keys, vals) )
        d['dist'] = d['displ']
        print('Units have been set to %s_%s' % (d['dist'], d['force']))        
        for k in d:
            print( 'For %s use %s' % (k, d[k]) )
        s.units = d
        
    def get_units(s, *args):
        # returs units for the passed argumnets
        
        res = []
        for a in args:
            mystr = ''
            if a in s.units:
                # check if the requested item is in the units dict
                mystr = ' ('+ s.units[a] + ')'
            elif a in _fieldtype:
                # return units on results dict item, ex: stress, strain etc
                ftype = _fieldtype[a]
                if ftype in s.units:
                    mystr = ' ('+ s.units[ftype] + ')'
            res.append(mystr)
        return res
        
    def set_time(s, time):
        # sets the time in the feamodel (preprocessor)
        s.time = time
        
    def get_item(s, item):
        # return item from part given string identifying item
        if item[0] == 'P':
            # get point
            items = s.points
            num = int(item[1:])
            res = [a for a in items if a.id == num]
            return res[0]
        elif item[0] == 'L':
            # get line
            items = s.lines
            num = int(item[1:])
            items = [a for a in items if a.id == num]
            return items[0]
        elif item[0] == 'A':
            # get area
            items = s.areas
            num = int(item[1:])
            items = [a for a in items if a.id == num]
            return items[0]
        else:
            print('Unknown item! Please pass the name of a point, line or area!')

    def MatlMaker(s, name):
        mat = Matl(name)
        s.matls.append(mat)
        return mat

    def PartMaker(s):
        p = Part(s)
        s.parts.append(p)
        return p

    def ModelMaker(s, parts, mtype):
        m = Model(s, parts, mtype)
        s.models.append(m)
        return m    
    
    def listify(s, items):
        # converts item into a list if it's not one
        if not isinstance(items, list):
            items = [items]
        return items
    
    def get_cname(s, items):
        # this returns a component name prefix, for labeling lists of items
        cname = ''
        if len(items) == 1:
            cname = items[0].get_name()
        else:
            cname = items[0].get_name()+'-'+items[-1].get_name()
        return cname

    def add_load(s, load, time):
        # add load to feamodel and part
        if time in s.loads:
            s.loads[time].append(load)
        else:
            s.loads[time] = [load]

    def set_gravity(s, grav, items):
        # applies gravity to the items
        items = s.listify(items)
        ctype = 'elements'        
        cname = s.get_cname(items)
        comp = Component(items, ctype, cname)
        
        if comp not in s.components:
            s.components.append(comp)
        
        ltype = 'gravity'
        time = s.time
        load = Load(ltype, comp, grav)
        
        # add load to feamodel
        s.add_load(load, time)

    def set_rpm(s, rpm, items):
        # applies rpm to the items
        items = s.listify(items)
        ctype = 'elements'
        cname = s.get_cname(items)
        comp = Component(items, ctype, cname)
        
        if comp not in s.components:
            s.components.append(comp)
        
        ltype = 'rpm'
        time = s.time
        load = Load(ltype, comp, rpm)
        
        # add load to feamodel
        s.add_load(load, time)

    def set_radps(s, rpm, items):
        # applies radps to the items
        items = s.listify(items)
        ctype = 'elements'
        cname = s.get_cname(items)
        comp = Component(items, ctype, cname)
        
        if comp not in s.components:
            s.components.append(comp)
        
        ltype = 'radps'
        time = s.time
        load = Load(ltype, comp, rpm)
        
        # add load to feamodel
        s.add_load(load, time)
        
    def set_fluid_press(s, str_items, rho, g, xo, po):
        # this sets pressure from water + atmosphere
        # p = po + rho*g*h, it assumes pressure increases in -x

        items = str_items
        # convert string into item(s)
        if isinstance(str_items, str):
            items = s.get_item(str_items)
        items = s.listify(items)

        ctype = 'faces'
        cname  = s.get_cname(items)
        comp = Component(items, ctype, cname)
        
        ltype = 'press_fluid'
        mult = rho*g
        load = Load_linear(ltype, comp, po, mult, xo)
        
        # add load to feamodel and part
        s.add_load(load, s.time)            

    def set_load(s, ltype, str_items, lval, ldir=None):
        # applies a load or pressure to lines
        # ltype = 'press' or 'force'
        # str_items is a string defining a line, point or side of the part
        #    or an item or a list of items
        # ldir = 'x' or 'y' or 'z'
        # positive pressure sign is compressing the surface
        # this is used to apply them to the child nodes/faces
        
        # make component if it doesn't exist
        
        # convert string into item(s)
        items = str_items
        # convert string into item(s)
        if isinstance(str_items, str):
            items = s.get_item(str_items)
        items = s.listify(items)
        
        ctype = 'nall'
        if ltype == 'press':
            ctype = 'faces'
        cname  = s.get_cname(items)
        comp = Component(items, ctype, cname)
        
        # check if component exists, and if so use it, need to write it here
        s.components.append(comp)

        if ltype == 'force':
            ltype = 'f'+ldir # for example fx

        # add load to feamodel
        load = Load(ltype, comp, lval)
        s.add_load(load, s.time)

    def set_constr(s, ltype, line_list, ldir, lval=0.0):
        # this sets constraints on lines
        # ltype = 'fix' or 'displ'
        # ldir = 'x' or 'y' or 'z'
        ctype = 'nall'

        # convert string into item(s)
        if isinstance(line_list, str):
            line_list = s.get_item(line_list)
        line_list = s.listify(line_list)

        cname  = s.get_cname(line_list)
        comp = Component(line_list, ctype, cname)

        # check if component exists, and if so use it, need to write it here
        s.components.append(comp)

        ltype = 'u'+ldir # for example ux
        load = Load(ltype, comp, lval)

        # add load to feamodel and part
        s.add_load(load, s.time)

    def set_eshape(s, eshape='quad', eorder=2):
        # sets eleemnt properties and thickness if applicable
        s.eshape = eshape # quad or tri
        s.eorder = eorder # 1 or 2

    def set_etype(s, items, etype='plstress', thick=None):
        # sets the element type on part areas, or a list of areas

        # convert string into item(s)
        if isinstance(items, str):
            items = s.get_item(items)

        # convert items into areas
        items = s.listify(items)
        cname = s.get_cname(items)
        if isinstance(items[0], Part):
            tmp = []
            for ind, part in enumerate(items):
                tmp += part.areas
            items = tmp

        # set the element types on the areas, ths is used to fix
        # elements when importing them from the inp file
        for area in items:
            area.set_etype(etype)

        # manually set all elments to the correct etype if they already exist
        for area in items:
            if hasattr(area, 'elements'):
                for e in area.elements:
                    e.set_etype(etype)

        # set a thickness component if needed
        if etype != 'axisym' and thick != None:
            
            # make component for nodal thickness
            ctype = 'nodes'
            comp = next((c for c in s.components if c.name == cname+'_nodes'), None)
            if comp == None:
                comp = Component(items, ctype, cname)
                s.components.append(comp)            

            # check existing thickness components and remove new area from them
            for time in s.loads:
                for load in s.loads[time]:
                    if load.ltype == 'thickness':
                        c = load.comp
                        areas = c.items
                        for newarea in items:
                            if neware in areas:
                                areas.remove(newarea)
                        # relabel modified component
                        c.name = s.get_cname(areas)+'_nodes'                        
                        
            ltype = 'thickness'
            time = 0.0
            load = Load(ltype, comp, thick)
            s.add_load(load, time)

    def set_matl(s, matl, items):
        # sets the matl on an item, part, or area + makes a component to apply it
        items = s.listify(items)
        cname = s.get_cname(items)
        if isinstance(items[0], Part):
            tmp = []
            for ind, part in enumerate(items):
                tmp += part.areas
            items = tmp

        # set the matls on the areas
        for area in items:
            pass
            #area.set_matl(matl)

        # make component
        ctype = 'elements'
        comp = next((c for c in s.components if c.name == cname+'_elements'), None)
        if comp == None:
            comp = Component(items, ctype, cname)
            s.components.append(comp)            

        # check existing matl components and remove new area from them
        for time in s.loads:
            for load in s.loads[time]:
                if load.ltype == 'matl':
                    c = load.comp
                    areas = c.items
                    for newarea in items:
                        if neware in areas:
                            areas.remove(newarea)
                    # relabel modified component
                    c.name = s.get_cname(areas)+'_elements'                        
        
        # manually set all elments to the correct matl if they already exist
        for area in items:
            if hasattr(area, 'elements'):
                for e in area.elements:
                    pass
        
        ltype = 'matl'
        time = 0.0
        load = Load(ltype, comp, matl)
        s.add_load(load, time)

    def mesh(s, fineness, mesher='cgx'):
        # this initiates meshing of all parts
        s.mesher = mesher
                
        if mesher == 'gmsh':
            s.mesh_gmsh(fineness)
        elif mesher == 'cgx':
            s.mesh_cgx(fineness)
    def mesh_gmsh(s, fine):
        # this meshes the model
        
        geo = []

        # write all points
        for pt in s.points:
            linestr = 'Point(%i) = {%f, %f, %f};' % (pt.id, pt.x, pt.y, 0.0)
            geo.append(linestr)
        
        # write all lines
        for line in s.lines:
            ln = line.id
            p1 = line.pt(0).id
            p2 = line.pt(1).id
            linestr = ''
            if isinstance(line, Arc):
                # line is arc
                pc = line.actr.id
                linestr = 'Circle(%i) = {%i, %i, %i};' % (ln, p1, pc, p2)
            else:
                # straight line
                linestr = 'Line(%i) = {%i,%i};' % (ln, p1, p2)
            geo.append(linestr)
            
            # set division if we have it
            if hasattr(line, 'ediv'):
                ndiv = line.ediv+1
                esize = line.length()/line.ediv
                if s.eshape == 'quad':
                    ndiv = line.ediv/2+1
                    esize = esize*2
                    # this is needed because quad recombine
                    # splits 1 element into 2
                linestr = 'Transfinite Line{%i} = %i;' % (ln, ndiv)
                print('LINE ELEMENT SIZE: %f, MAKES %i ELEMENTS' % (line.length()/line.ediv, line.ediv))
                geo.append(linestr)
                geo.append('Characteristic Length {%i,%i} = %f;' % (p1, p2, esize))
            

        # write all areas
        for area in s.areas:
            if area.closed:
                aid = area.id
                aname = area.get_name()
                linestr = 'Line Loop(%i) = ' % (aid)
                line_ids = []
                for line in area.lines:
                    char = ''
                    if line.sign == -1:
                        char = '-'
                    line_ids.append(char+str(line.id))
                linestr = linestr + '{'+ ','.join(line_ids)+'};'
                geo.append(linestr)
                geo.append('Plane Surface(%i) = {%i};' % (aid, aid))
                geo.append("Physical Surface('%s') = {%i};" % (aname, aid))

        # write part area components
        for part in s.parts:
            # make components for each part
            line = "Physical Surface('%s') = " % (part.get_name())
            area_ids = []
            for area in part.areas:
                if area.closed:
                    area_ids.append(str(area.id))
            line = line + '{' + ','.join(area_ids) + '};'
            geo.append(line)
                    
        # write all line componenets so we can get nodes out
        for L in s.lines:
            line = "Physical Line('%s') = {%i};" % (L.get_name(), L.id)
            geo.append(line)
        
        # write node componenets
        # node list is not produced by gmsh
        for pt in s.points:
            linestr = "Physical Point('%s') = {%i};" % (pt.get_name(), pt.id)
            geo.append(linestr)                    

        # set the meshing options
        geo.append('Mesh.CharacteristicLengthFactor = '+str(fine)+'; //mesh fineness')
        geo.append('Mesh.RecombinationAlgorithm = 1; //blossom')
        #geo.append('Mesh.Lloyd = 1; //smoothing algorithm')

        if s.eshape == 'quad':
            geo.append('Mesh.RecombineAll = 1; //turns on quads')
            geo.append('Mesh.SubdivisionAlgorithm = 1; // quadrangles only')
            #geo.append('Mesh.RecombinationAlgorithm = 1; //turns on blossom needed for quad')

        eo = s.eorder
        geo.append('Mesh.CharacteristicLengthExtendFromBoundary = 1;')
        geo.append('Mesh.CharacteristicLengthMin = 0;')
        geo.append('Mesh.CharacteristicLengthMax = 1e+022;')
        # use this so small circles are meshed finely
        geo.append('Mesh.CharacteristicLengthFromCurvature = 1;')
        geo.append('Mesh.CharacteristicLengthFromPoints = 1;')
        #geo.append('Mesh.Algorithm = 2; //delauny') #okay for quads
        geo.append('Mesh.Algorithm = 8; //delquad = delauny for quads')
        geo.append('Mesh.ElementOrder = '+str(eo)+'; //linear or second set here')
        if eo == 2:
            geo.append('Mesh.SecondOrderIncomplete=1; //req for 2nd ord, no face node')
        geo.append('Mesh.SaveGroupsOfNodes = 1; // save node groups')
                
        # write geo file to the local directory
        fname = s.fname+'.geo'
        fout = s.fname+'.inp'
        f = open(fname,'w')
        for line in geo:
            #print (line)
            f.write(line+'\n')
        f.close()
        print ('File: '+ fname + ' was written')
        
        # run file in bg mode, -2 is 2d mesh
        #p = subprocess.Popen("%s %s -2 -order %i -o %s" % (_gmsh, fname, eo, fout))
        p = subprocess.Popen("%s %s -2 -o %s" % (_gmsh, fname, fout))
        p.wait()
        print ('File: '+ fout + ' was written')
        print('Meshing done!')
        
        # fix element type
        # define the element types, and tells the program to mesh
        estr = s.eshape+str(s.eorder)+s.parts[0].areas[0].etype
        ccxe = _ccx_elements[estr]
        #print('PREDICTED CCX ETYPE IS: ',ccxe)
        
        # change the part elements to the right type
        inp = []
        f = open(fout,'r')
        for line in f:
            if '*Element, type=' in line:
                etype = line[15:19]
                if 'CPS' in etype:
                    # fix element type here, it will be CPS3,4,6 or 8
                    # replace PS with ES or AX
                    estr = ccxe[1:3]
                    line = line[:16] + estr + line[18:]
            line = line.strip()
            inp.append(line)
        f.close()
        # write fixed file
        f = open(fout,'w')
        for line in inp:
            #print (line)
            f.write(line+'\n')
        f.close()

        # write gmsh msh file
        p = subprocess.Popen("%s %s -2 -o %s" % (_gmsh, fname, s.fname+'.msh'))
        #p = subprocess.Popen("%s %s -2 -o %s" % (_gmsh, fname, s.fname+'.jpg'))
        print ('File: '+ s.fname+ '.msh was written')

        # open the output file if it is a gmsh msh file
        if fout[-3:] == 'msh':
            # start gmsh if fout file is gmesh mesh
            p = subprocess.Popen("%s %s" % (_gmsh, fout))
            p.wait()

        # read in the calculix mesh
        s.read_inp(s.fname+'.inp')

    def read_inp(s, fname):
        # this function reads in the calculix input file
        # stores all nodes, elements, and faces
        # and assigns all element and node sets correctly

        f = open(fname,'r')
        mode = None
        set_name = None
        set_type = None
        
        items = [] # holder for nodes or elements in nsets or esets
        N = ID_List() # store nodes        
        E = ID_List() # store elements, allows renumbering before putting int model
        F = [] # store faces
        sets = {'E':{},'N':{}} # store sets
        
        # read in input file
        for line in f:
            if line[0] != '*':
                if mode == 'nmake':
                    L = line.split(',')
                    L = [a.strip() for a in L]
                    (nnum, x, y, z) = (int(L[0]), float(L[1]), float(L[2]), float(L[3]))
                    node = Node(nnum, x, y, z)
                    N.append(node)
                elif mode == 'emake':
                    L = line.split(',')
                    L = [int(a.strip()) for a in L]
                    enum = L[0]
                    nlist = [N.idget(a) for a in L[1:]]
                    e = Element(enum, etype, nlist)
                    faces = e.faces()
                    E.append(e)
                    F += faces
                    sets[set_type][set_name].append(e)
                elif mode == 'set':
                    L = line.split(',')
                    L = [a.strip() for a in L]
                    L = [int(a) for a in L if a != '']
                    items = []
                    if set_type == 'E':
                        items = [E.idget(a) for a in L]
                    elif set_type == 'N':
                        items = [N.idget(a) for a in L]
                    if items == [None]*len(items):
                        pass # the elements were not found
                    else:
                        sets[set_type][set_name] += items
                    
            # mode setting
            if '*Node' in line or '*NODE' in line:
                mode = 'nmake'
            elif '*Element' in line or '*ELEMENT' in line:
                L = line.split(',') # split it based on commas
                e = L[1].split('=')
                etype = e[1]
                
                # exclude T elements made in gmsh
                if etype[0] != 'T':
                    e = L[2].split('=')
                    set_name = e[1].strip()
                    set_type = 'E'
                    sets[set_type][set_name] = []
                    mode = 'emake'
                else:
                    mode = None
            elif '*ELSET' in line:                    
                L = line.split(',')
                e = L[1].split('=')
                set_name = e[1].strip()
                set_type = 'E'
                sets[set_type][set_name] = []
                mode = 'set'
            elif '*NSET' in line:        
                L = line.split(',')
                e = L[1].split('=')
                set_name = e[1].strip()
                set_type = 'N'
                sets[set_type][set_name] = []
                mode = 'set'
        f.close()        
        
        # loop through sets and remove empty sets
        # store sets to delete
        todel = []
        for (set_type, set_dict) in sets.items():
            for (set_name, item_list) in set_dict.items():
                if item_list == []:
                    todel.append( {'set_type':set_type,'set_name':set_name} )
        # delete the empty sets
        for d in todel:
            (set_type, set_name) = (d['set_type'],d['set_name'])
            del sets[set_type][set_name]
            #print('Empty set type:%s name:%s deleted' % (set_type, set_name))

        # loop through faces and set externa flag to true
        external = []
        for face in F:
            num = F.count(face)
            if num == 1:
                external.append(face)
        for f in external:
            # set external flag to true
            f.set_ext()

        # this resets the min element to number 1
        if E.get_minid() > 1:
            E.set_minid(1)

        #-----------------------------------
        # Node and element assignment back onto parts, areas, lines, points
        #-----------------------------------
        # assign elements + nodes to parts
        s.elements = E
        s.nodes = N
        s.faces = F
        
        for part in s.parts:
            # assign part element and node sets
            pname = part.get_name()
            part.elements = sets['E'][pname]
            part.nodes = sets['N'][pname]
            
            # assign all nodes and elements to areas 
            for area in s.areas:
                aname = area.get_name()
                # change the element to the right type based on python type
                for ind in range(len(sets['E'][aname])):
                    sets['E'][aname][ind].set_etype( area.etype )
                area.elements = sets['E'][aname]
                area.nodes = sets['N'][aname]
            
            # assign the child nodes to points
            pts = part.get_points()
            for pt in pts:
                ndist = []
                for n in part.nodes:
                    p_tmp = Point(n.x, n.y)
                    p_tmp = pt - p_tmp
                    dist = p_tmp.length()
                    ndist.append( {'dist':dist,'node':n} )
                # sort the list by dist, sorts low to high
                ndist = sorted(ndist, key=lambda k: k['dist'])
                pt.nall = ndist[0]['node']
                print('Point %s = node %s' % (pt, pt.nall))
            
            # assign the nall and n1 and faces to lines
            lines = part.get_lines()
            for line in lines:
                lname = line.get_name()
                nall = sets['N'][lname]
                n1 = [n for n in nall if n.order == 1]
                faces = []
                for face in external:
                    if set(face.nodes).issubset(set(nall)):
                        faces.append(face)
                line.nall = nall
                line.n1 = n1
                line.faces = faces
                
        print('Done reading Calculix/Abaqus .inp file')


    def mesh_cgx(s, fine):
        # meshes the parts using calculix preprocessor
        fbd = []
        comps = []
        cfiles = []
        
        num = 1.0/fine
        emult = int(round(num)) # this converts fineness to mesh multiplier
        
        # write all points
        for pt in s.points:
            linestr = 'pnt %s %f %f %f' % (pt.get_name(), pt.x, pt.y, 0.0)
            fbd.append(linestr)
            # gmsh can't make node componenets so don't do it in cgx
            #L = 'seta %s p %s' % (pt.get_name(), pt.get_name())
            #comps.append(L)            
        
        # write all lines
        for line in s.lines:
            ln = line.get_name()
            p1 = line.pt(0).get_name()
            p2 = line.pt(1).get_name()
            linestr = ''
            if isinstance(line, Arc):
                # line is arc
                pc = line.actr.get_name()
                linestr = 'line %s %s %s %s' % (ln, p1, p2, pc)
            else:
                # straight line
                linestr = 'line %s %s %s' % (ln, p1, p2)
            # set division if we have it
            if hasattr(line, 'ediv'):
                for part in s.parts:
                    if line in part.get_lines():
                        break
                ndiv = part.eorder*line.ediv 
                linestr += ' '+str(int(ndiv))
            fbd.append(linestr)
            L = 'seta %s l %s' % (ln, ln)
            comps.append(L)
            cfiles.append(ln)
        
        # write all areas
        for area in s.areas:
            if area.closed:
                linestr = 'gsur '+area.get_name()+' + BLEND '
                line_ids = []
                for line in area.lines:
                    char = '+'
                    if line.sign == -1:
                        char = '-'
                    line_ids.append(char+' '+line.get_name())
                linestr = linestr + ' '.join(line_ids)
                fbd.append(linestr)
                # add area component, nodes + elements
                L = 'seta %s s %s' % (area.get_name(), area.get_name())
                comps.append(L)
                cfiles.append(area.get_name())
                
        
        # write part area components
        for part in s.parts:
            # make components for each part, syntax
            # seta P0 s s0 s1
            line = 'seta %s s ' % (part.get_name(),)
            cfiles.append(part.get_name())
            area_ids = []
            for area in part.areas:
                if area.closed:
                    area_ids.append(area.get_name())
            line = line + ' '.join(area_ids)
            fbd.append(line)
            #lines that mesh the part
            part_comp = part.get_name()
            
            # define the element types, and tells the program to mesh
            estr = part.eshape+str(part.eorder)+part.etype
            etype = _cgx_elements[estr]
            fbd.append('elty '+part_comp+' '+etype)
            fbd.append('div all mult '+str(emult))
            #fbd.append('div all auto 2. 10. 0.5')
            fbd.append('mesh '+part_comp)
            # this is needed to select all nodes of the elements

        # save mesh file
        fbd.append('send all abq')        
        
        # add line and area components
        fbd += comps
        
        # save component node and element sets
        for comp in cfiles:
            # select nodes under
            fbd.append('comp %s do' % (comp,))
            fbd.append('send %s abq names' % (comp,))
        
        # this orients the view correctly
        # this is the same as switching to the z+ orientation
        # y us axial, x is radial
        fbd.append('rot z')
        fbd.append('rot c -90')
        fbd.append('plot e all')
        
        # write fbd file to the local directory
        fname = s.fname+'.fbd'
        f = open(fname,'w')
        for line in fbd:
            #print (line)
            f.write(line+'\n')
        f.close()
        print ('File: '+ fname + ' was written')
        
        # run file in bg mode
        p = subprocess.Popen("%s -bg %s" % (_cgx, fname))
        p.wait()
        print('Meshing done!')
        
        # assemble the output files into a ccx input file
        inp = []
        files = ['all.msh']
        files += [f+'.nam' for f in cfiles]
        for fname in files:
            f = open(fname,'r')
            for line in f:
                # cgx adds E and N prfixes on sets after =, get rid of these
                if '=' in line and fname != 'all.msh':
                    L = line.split('=')
                    line = L[0] + '=' + L[1][1:]
                    inp.append(line.strip())
                else:
                    inp.append(line.strip())
            f.close()
            
            # delete file
            os.remove(fname)
        
        # write out inp file
        fname = s.fname+'.inp'
        f = open(fname,'w')
        for line in inp:
            #print (line)
            f.write(line+'\n')
        f.close()
        print ('File: '+ fname + ' was written')

        # read in the calculix mesh
        s.read_inp(fname)

class Matl(Idobj):
    # stores a material
    def __init__(s, name):
        s.name = name
        Idobj.__init__(s)        
    def set_mech_props(s, density, youngs, pratio):
        s.density = density
        s.youngs = youngs
        s.pratio = pratio
    def set_therm_props(s, conductivity, spec_heat):
        s.conductivity = conductivity
        s.spec_heat = spec_heat
    def set_therm_expan(s, alphas, temps, tzero):
        s.thermal_exp = zip(alphas, temps)
        s.tzero = tzero
    def ccx(s):
        # this writes out a material in ccx form
        res = []
        res.append('*MATERIAL,NAME='+s.name)
        if hasattr(s,'youngs'):
            res.append('*ELASTIC')
            res.append(str(s.youngs)+','+str(s.pratio))
        if hasattr(s, 'density'):
            res.append('*DENSITY')
            res.append(str(s.density))
        if hasattr(s, 'conductivity'):
            res.append('*CONDUCTIVITY')
            res.append(str(s.conductivity))
        if hasattr(s, 'spec_heat'):
            res.append('*SPECIFIC HEAT')
            res.append(str(s.spec_heat))
        if hasattr(s, 'thermal_exp'):
            res.append('*EXPANSION,ZERO='+str(s.tzero))
            for pair in s.thermal_exp:
                res.append(str(pair[0])+' '+str(pair[1]))
        return res