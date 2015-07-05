"""This module stores geometry classes, which are used to make parts.
"""
from math import atan2, pi, cos, sin, radians, degrees
from matplotlib.patches import Arc as AArc
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from numpy.linalg import det
from numpy import array
from numpy.linalg import solve as linsolve

from . import base_classes

#accuracy for small numbers in math below
ACC = .00001 # original

LDIVS = {}

def get_text_hv(angle):
    """Returns (ha, va) text horizontal and vertical alignment for line label.

    This makes the line label text inside its area, assuming areas closed CW.

    Args:
        angle (float): line normal vector in degrees, between 180 and -180
    """
    horiz, vert = ('', '')
    if abs(angle) < 67.5:
        horiz = 'right'
    elif 67.5 <= abs(angle) < 112.5:
        horiz = 'center'
    elif 112.5 <= abs(angle) <= 180:
        horiz = 'left'
    if abs(angle) < 22.5 or abs(angle) > 157.5:
        vert = 'center'
    else:
        if angle > 0:
            vert = 'top'
        else:
            vert = 'bottom'
    return (horiz, vert)


class Point(base_classes.Idobj):
    """Makes a point or vector.

    Args:
        x (float): vertical coordinate of point
        y (float): horizontal coordinate of point
        z (float): in-page coordinate of point, defaults to zero

    Attributes:
        x (float): vertical coordinate of point
        y (float): horizontal coordinate of point
        z (float): in-page coordinate of point, defaults to zero
            nodes (Node): a list of nodes in a mesh at this point.
            List length is 1.
        nodes (list): list of nodes under point
        lines (list): list of lines which use this point
        arc_center (bool): True if this point is an arc center
    """

    def __init__(self, x, y, z=0):
        self.x = x
        self.y = y
        self.z = z
        self.nodes = []
        self.lines = []
        self.arc_center = False
        self.on_part = True
        base_classes.Idobj.__init__(self)

    def get_name(self):
        """Returns item name."""
        return 'P'+str(self.id)

    def __eq__(self, other):
        """Checks x and y equality to other point.

        Args:
            other (Point): passed point
        """
        if isinstance(other, Point):
            if self.x == other.x and self.y == other.y:
                return True
            return False
        else:
            return False

    def __add__(self, other):
        """Returns new point = self + other.

        Args:
            other (Point): passed point

        Returns:
            Point = self + other (new point, does not modify self)
        """
        return Point(self.x+other.x, self.y+other.y, self.z+other.z)

    def __sub__(self, other):
        """Returns new point = self - other.

        Args:
            other (Point): passed point

        Returns:
            Point = self - other (new point, does not modify self)
        """
        return Point(self.x-other.x, self.y-other.y, self.z-other.z)

    def __hash__(self):
        """Set the hash of the item to the id number.

        This allows one to make sets of these items.
        """
        return self.id

    def __mul__(self, factor):
        """Returns new point = self*factor.

        Args:
            factor (float): factor to multiply the Point by

        Returns:
            Point = self*factor (new point, does not modify self)
        """
        return Point(self.x*factor, self.y*factor, self.z*factor)

    def __truediv__(self, other):
        """Returns new point = self/other.

        Args:
            other (float or Point): item we're dividing by

        Returns:
            Point or factor. Does not modify self.
        """
        if isinstance(other, float):
            return Point(self.x/other, self.y/other, self.z/other)
        elif isinstance(other, Point):
            res = []
            if other.x != 0:
                res.append(self.x/other.x)
            if other.y != 0:
                res.append(self.y/other.y)
            if other.z != 0:
                res.append(self.z/other.z)
            factor = sum(res)/len(res)
            return factor

    def save_line(self, line):
        """Saves the line or arc in the point."""
        if line not in self.lines:
            self.lines.append(line)
        if self.arc_center == True:
            if len(self.lines) > 1:
                for line in self.lines:
                    if isinstance(line, Line):
                        self.arc_center = False
                        break

    def length(self):
        """Returns the length of this point or vector."""
        res = (self.x**2 + self.y**2)**0.5
        return res

    def label(self, axis):
        """Labels the point on a Matplotlib axis.

        Args:
            axis (Matplotlib Axis): Matplotlib Axis
        """
        axis.annotate(self.get_name(), (self.y, self.x))

    def plot(self, axis, label=True):
        """Plots the point on the passed matplotlib axis.

        Args:
            axis (Matplotlib axis): plate to plot the point
            pnum (bool): True turns on point labeling
        """
        axis.scatter(self.y, self.x)
        if label:
            self.label(axis)

    def make_unit(self):
        """Modifies self vector to make it a unit vector."""
        lenval = self.length()
        self.x = self.x*1/lenval
        self.y = self.y*1/lenval

    def ang_rad(self):
        """Returns angle in radians, assume point is vector from 0,0."""
        vert = self.x
        horiz = self.y
        rad_angle = atan2(vert, horiz)
        return rad_angle

    def ang_deg(self):
        """Returns angle in degrees, assume point is vector from 0,0."""
        rad_angle = self.ang_rad()
        deg = rad_angle * 180.0 / pi
        return deg

    def axrad(self):
        """Returns a tuple of (ax, rad) for plotting."""
        return (self.y, self.x)

    def radax(self):
        """Returns a tuple of (rad, ax) for checking if inside."""
        return (self.x, self.y)

    def rot_ccw_deg(self, ang):
        """Rotates the current vector by ccw degrees about 0,0. Returns self.

        Args:
            ang (float): angle or rotation in degrees, counter-clockwise (ccw)
        """
        ang = radians(ang)
        axial = self.y*cos(ang) - self.x*sin(ang)
        rad = self.y*sin(ang) + self.x*cos(ang)
        axial = round(axial, 5)
        rad = round(rad, 5)
        (self.x, self.y) = (rad, axial)
        return self

    def ang_bet_rad(self, other):
        """Returns angle between self and other vector in radians.

        Assumes vectors start at 0,0

        Args:
            other (Point): passed vector

        Returns:
            angle (float): radians between self vector and other vector
        """
        avect = Point(self.y*other.x - self.x*other.y, self.y*other.y + self.x*other.x)
        ang = avect.ang_rad()
        return ang

    def ang_bet_deg(self, other):
        """Returns angle between self and other vector in degrees.

        Assumes vectors start at 0,0

        Args:
            other (Point): passed vector

        Returns:
            angle (float): degrees between self vector and other vector
        """
        avect = Point(self.y*other.x - self.x*other.y, self.y*other.y + self.x*other.x)
        ang = avect.ang_deg()
        return ang

    def __str__(self):
        """Returns string listing object type, id number, and coordinates"""
        val = 'Point %s, (x, y)=(%.3f, %.3f)' % (self.get_name(), self.x, self.y)
        return val


class Line(base_classes.Idobj):
    """Stores a line from points start_pt to end_pt.

    Args:
        start_pt (Point): first point
        end_pt (Point): second point

    Attributes:
        points (list of Point) : list of the line's points [start_pt, end_pt]
        ediv (None or float): number of elements on the line
        midpt (Point): a mid-point between start_pt and end_pt
        nodes (list of Node): a list of meshed nodes on the line
        edge (bool): True if line is a non-shared edge in an area
        signlines (list): list of SignLine that use this Line
    """

    def __init__(self, start_pt, end_pt):
        self.points = [start_pt, end_pt]
        self.midpt = self.mid()
        self.ediv = None
        self.nodes = []
        self.edge = True
        self.signlines = []
        base_classes.Idobj.__init__(self)

    def reverse(self):
        """Reverses self."""
        self.points.reverse()

    @property
    def allpoints(self):
        """Returns all line defining points."""
        return self.points

    @property
    def areas(self):
        """Returns all areas that signlines use."""
        areas = set()
        for sline in self.signlines:
            areas.add(sline.parent)
        return list(areas)

    def __eq__(self, other):
        """Returns bool telling if self == other.

        Returns True if the lines are the same type, and contain the same
        list of points.
        """
        if isinstance(other, self.__class__):
            mypts = self.points
            otherpts = other.points
            same = True
            for (mypt, otherpt) in zip(mypts, otherpts):
                if mypt != otherpt:
                    same = False
                    break
            return same
        else:
            return False

    def __hash__(self):
        """Set the hash of the item to the id number.

        This allows one to make sets of these items.
        """
        return self.id

    def add_signline(self, signline):
        """Adds a signline to self.signlines."""
        if signline not in self.signlines:
            self.signlines.append(signline)

    def get_name(self):
        """Returns the line name."""
        return 'L'+str(self.id)

    def length(self):
        """Returns the length of the line."""
        return (self.pt(0)-self.pt(1)).length()

    def save_to_points(self):
        """This stores this line in the line's child points."""
        for point in self.allpoints:
            point.save_line(self)

    def set_ediv(self, ediv):
        """Sets the number of element divisions on the line when meshing.

        Args:
            ediv (int): number of required elements on this line
        """
        self.ediv = ediv

    def signed_copy(self, sign):
        """Returns a SignLine copy of this Line with the passed sign."""
        return SignLine(self, sign)

    def pt(self, ind):
        """Returns the start or end point of the line.

        Args:
            ind (int): index of the point we want, 0=start, 1=end

        Returns:
            Point: requested Point
        """
        return self.points[ind]

    def set_pt(self, ind, point):
        """Update the line's point at ind to the new point.

        Args:
            ind (int): index of the point we're updating, 0=start, 1=end
            point (Point): new point we are assigning to Line
        """
        # remove this line from the old point
        self.points[ind].lines.remove(self)
        # update the point
        self.points[ind] = point
        self.midpt = self.mid()
        self.save_to_points()
        other_point = self.pt(int(not ind))
        if point == other_point:
            # we are making a zero length line
            other_point.lines.remove(self)
            if self.id != -1:
                if len(self.signlines) > 0:
                    area = self.signlines[0].lineloop.parent
                    area.part.fea.lines.remove(self)
            for sline in self.signlines:
                sline.lineloop.remove(sline)

    def mid(self):
        """Calculates and returns a midpoint between start and end points."""
        point = (self.pt(0) + self.pt(1))*0.5
        return point

    def touches(self, other):
        """Checks if a line is connected to a passed line.

        Args:
            other (Line or Arc): passed other line

        Returns:
            True: beginning or end point connect to other line
            False: lines are not connected
        """
        [start_pt, end_pt] = [self.pt(0), self.pt(1)]
        if start_pt in other.points or end_pt in other.points:
            return True
        else:
            return False

    def offset(self, dist):
        """Returns a new line offset from the dirrent line by value dist.

        Right hand rule is followed assuming x is the origina line.
        Z points out of page.
        Y is the direction of line offset.

        Args:
            dist (float): distance to offset the new line by

        Returns:
            tmpline (Line): new Line which is offset from self by dist
        """
        lvect = self.get_perp_vec()
        lvect.make_unit()
        # scale it by the new amount
        lvect = lvect*dist
        # add the vector onto the defining points
        start_pt = self.pt(0) + lvect
        end_pt = self.pt(1) + lvect
        tmpline = Line(start_pt, end_pt)
        return tmpline

    def get_abc(self):
        """ Returns a list of abc terms for a line definition.

        ax + by + c = 0

        Returns:
            [a, b, c]: list of the above line terms
        """
        delta = self.pt(1)-self.pt(0)
        dx_ = delta.x
        dy_ = delta.y
        (a__, b__, c__) = (0.0, 0.0, 0.0)
        if dy_ == 0:
            # vertical radial line
            # y = 5 --> 0 = -y + 5 --> 0 = 0x -y + 5
            (b__, c__) = (-1.0, self.pt(0).y)
        elif dx_ == 0:
            # horizontal axial line
            # x = 5 --> 0 = -x + 5 --> 0 = -x +0y + 5
            (a__, c__) = (-1.0, self.pt(0).x)
        else:
            # y = ax + c --> 0 = ax - y + c
            slope = dy_/dx_
            offset = self.pt(0).y - slope*self.pt(0).x
            (a__, b__, c__) = (slope, -1.0, offset)
        return [a__, b__, c__]

    def get_perp_vec(self, point=None):
        """ Returns vector perpendicular to current line.

        Vector is created by rotating the current line ccw 90 degrees
        around the start point.

        Args:
            point (Point): defaults to None, unused for Line, used for Arc

        Returns:
            lvect (Point): vector that is perpendiculat to the current line
        """
        lvect = self.pt(1) - self.pt(0)
        lvect.rot_ccw_deg(90)
        return lvect

    def get_tan_vec(self, point):
        """ Returns vector tangent to the current line.

        Vector points from the passed point to a unit location further away.

        Args:
            point (Point): start point

        Returns:
            lvect (Point): vector that is tangent to the current line
        """
        end = self.pt(0)
        if point == end:
            end = self.pt(1)
        lvect = end - point
        lvect.make_unit()
        return lvect

    def arc_tang_intersection(self, point, mag):
        """Returns an intersection point on this line of an arc at center point

        Args:
            point (Point): arc center point
            mag (float): passed radius distance from the arc center to the line

        Returns:
            newpt (Point or None):
                Point: the intersection point
                None: returns None if no intersection point exists
        """
        vect = self.get_perp_vec()
        vect.make_unit()
        pt2 = point + vect*(-2*mag)

        tmpline = Line(point, pt2)
        newpt = tmpline.intersects(self)
        if type(newpt) == type(None):
            print('Intersection failed!')
        return newpt

    def coincident(self, point):
        """Checks to see if point is on the line.

        Args:
            point (Point): input point to check

        Returns:
            bool: True if point on line or False if not
        """
        # check if point is almost equal to end points
        dist = point - self.pt(0), point - self.pt(1)
        dist = [item.length() for item in dist]
        if min(dist) < ACC:
            return True
        # check if point is on the lines
        [a__, b__, c__] = self.get_abc()
        remainder = abs(a__*point.x + b__*point.y + c__)
        # print('  coincident: %s %s' % (self.get_name(), point))
        # print('  remainder: %.6f' % remainder)
        if remainder < ACC:
            # point is on the line equation, but is it between end pts
            xvals = [apoint.x for apoint in self.points]
            yvals = [apoint.y for apoint in self.points]
            xin = min(xvals)-ACC <= point.x <= max(xvals)+ACC
            yin = min(yvals)-ACC <= point.y <= max(yvals)+ACC
            if xin and yin:
                # we're on the line
                return True
            else:
                # we're off of the line
                # print('  outside of bounds!')
                return False
        else:
            # point is not on line equation
            return False

    def intersects(self, other):
        """Checks if other line intersects this line.

        Args:
            other (Line, Arc, Signline, SignArc): other line to check

        Returns:
            Point or None:
                Point: intersection point
                None: None if lines don't intersect
        """
        # other is cutline when checking
        if isinstance(other, Line) or isinstance(other, SignLine):
            terms_1 = self.get_abc()
            terms_2 = other.get_abc()
            p_1, p_2 = self.pt(0), self.pt(1)
            p_3, p_4 = other.pt(0), other.pt(1)
            dist = [p_1 - p_3, p_1 - p_4, p_2 - p_3, p_2 - p_4]
            vals = [p_1, p_1, p_2, p_2]
            dist = [item.length() for item in dist]
            # check each of self's point to see if it == other's points
            for (distance, val) in zip(dist, vals):
                if distance < ACC:
                    return Point(val.x, val.y)
            if terms_1[:2] == terms_2[:2]:
                # lines are parallel and will not have an intersection point
                #print('Line-Line intersection = None: same slope')
                return None
            denom = (p_1.x-p_2.x)*(p_3.y-p_4.y) - (p_1.y-p_2.y)*(p_3.x-p_4.x)
            if abs(denom) < ACC:
                # lines are nearly identical slope
                #print('Line-Line intersection = None because of denom')
                return None
            else:
                # intersection calculation
                a_array = array([[terms_1[0], terms_1[1]],
                                 [terms_2[0], terms_2[1]]])
                b_array = array([-terms_1[2], -terms_2[2]])
                result = linsolve(a_array, b_array)
                point = Point(result[0], result[1])
                # check against existing points on self
                dist = [point - p_1, point - p_2]
                vals = [p_1, p_2]
                dist = [item.length() for item in dist]
                for (distance, loc) in zip(dist, vals):
                    if distance < ACC:
                        point = Point(loc.x, loc.y)
                        break
                on_self = self.coincident(point)
                on_other = other.coincident(point)
                if on_self and on_other:
                    """
                    s_name = self.get_name()
                    o_name = other.get_name()
                    print('X Checking intersection of %s and %s' % (s_name, o_name))
                    print('  %s' % point)
                    print('  (on_self, on_other): (%s, %s)' % (on_self, on_other))
                    """
                    return point
                else:
                    """
                    print('Line-Line intersection = None point not on lines')
                    print(' %s' % point)
                    print(' on_self: %s' % on_self)
                    print(' on_other: %s' % on_other)
                    """
                    return None
        elif isinstance(other, Arc) or isinstance(other, SignArc):
            # arc line intersection
            return other.instersects(self)

    def __str__(self):
        """Returns string listing object type, name, and points"""
        start = self.pt(0)
        end = self.pt(1)
        val = 'Line %s, start: %s end: %s' % (self.get_name(), start, end)
        return val

class SignLine(Line, base_classes.Idobj):
    """Makes a signed line.

    Attributes:
        line (Line): parent line
        sign (int): 1 = positive, or -1 = negative
        lineloop (None or LineLoop): parent LineLoop, default is None
        faces (list): list of element faces
        n1 (list): list of first order nodes on line
        nodes (list): nodes on the line
    """
    def __init__(self, parent_line, sign):
        self.line = parent_line
        self.sign = sign
        self.lineloop = None
        self.faces = []
        self.n1 = []
        base_classes.Idobj.__init__(self)

    def pt(self, index):
        """Returns the point at the start or end of the line.

        Args:
            index (int): 0 = start point, 1 = end point
        """
        if self.sign == -1:
            return self.line.points[not index]
        else:
            return self.line.points[index]

    def set_pt(self, index, point):
        """Sets the start or end point of the line to a new point.

        Args:
            index (int): 0 for start point, 1 for end point
            point (Point): the new point we are updating to
        """
        self.line.set_pt(index, point)

    def set_ediv(self, ediv):
        """Applies the element divisions onto the parent line."""
        self.line.set_ediv(ediv)

    def set_lineloop(self, lineloop):
        """Sets the parent LineLoop"""
        # this is needed to cascade set ediv up to FEA model and down onto
        # other lines with this ediv
        self.lineloop = lineloop

    def get_name(self):
        """Returns the SignLine name."""
        if self.sign == -1:
            return '-L'+str(self.line.id)
        else:
            return 'L'+str(self.line.id)

    def label(self, ax):
        """Labels the line on the matplotlib ax axis.

        Args:
            ax (Matplotlib axis): Matplotlib axis
        """
        # find angle of perpendiclar line so we can set text alignment
        lvect = self.get_perp_vec()
        ang = lvect.ang_deg()
        horiz, vert = get_text_hv(ang)
        axial = self.midpt.y
        radial = self.midpt.x
        ax.text(axial, radial, self.get_name(), ha=horiz, va=vert)

    def plot(self, ax, label=True):
        """Draws and labels the line onto the passed Matplotlib axis.

        Args:
            ax (Matlotlib axis): Matplotlib axis
            label (bool): if True label the line
        """
        lax = [pt.y for pt in self.points]
        lrad = [pt.x for pt in self.points]
        # draw line
        ax.plot(lax, lrad)
        if label:
            self.label(ax)

    def signed_copy(self, sign):
        """Returns a SignLine copy of this line."""
        realsign = self.sign*sign
        return SignLine(self.line, realsign)

    @property
    def nodes(self):
        """Gets the line nodes."""
        return self.line.nodes

    @nodes.setter
    def nodes(self, nlist):
        """Sets the line nodes."""
        self.line.nodes = nlist

    @property
    def midpt(self):
        """Returns the mid point."""
        return self.line.midpt

    @property
    def edge(self):
        """Returns a bool saying if this is an edge line."""
        return self.line.edge

    @edge.setter
    def edge(self, isedge):
        """Sets the line nodes."""
        self.line.edge = isedge

    @property
    def points(self):
        """Returns a list of the line's points [pstart, pend]"""
        pstart = self.pt(0)
        pend = self.pt(1)
        return [pstart, pend]

    def reverse(self):
        """Reverses self.line"""
        self.line.reverse()


class Arc(base_classes.Idobj):
    """Makes an arc from points start_pt to end_pt about center actr.

    Arcs should be 90 degrees maximum.

    Arguments:
        start_pt (Point): first point
        end_pt (Point): second point
        actr (Point): arc center

    Attributes:
        points (list of Point) : list of the line's points [start_pt, end_pt]
        allpoints (list of Point): [start_pt, end_pt, actr]
        actr (Point): arc center
        radius (float): radius of the arc
        concavity (str):

            - 'concave': concave arc
            - 'convex': convex arc

        midpt (Point): a mid-point between p1 and p2
        ediv (None or int): number of elements on arc
        nodes (list of Node): a list of meshed nodes on the line
        edge (bool): True if arc is a non-shared edge in an area
        signlines (list): list of SignArc that use this Arc
    """

    def __init__(self, start_pt, end_pt, actr):
        self.points = [start_pt, end_pt]
        self.actr = actr
        self.actr.arc_center = True
        self.radius = (start_pt - actr).length()
        self.concavity = self.get_concavity()
        self.midpt = self.mid()
        self.ediv = None
        self.nodes = []
        self.edge = True
        self.signlines = []
        base_classes.Idobj.__init__(self)

    def reverse(self):
        """Reverses self."""
        self.points.reverse()
        self.concavity = self.get_concavity()

    @property
    def allpoints(self):
        """Returns all points including arc center"""
        return self.points + [self.actr]

    @property
    def areas(self):
        """Returns all areas that signlines use."""
        areas = set()
        for sline in self.signlines:
            areas.add(sline.parent)
        return list(areas)

    def __eq__(self, other):
        """Returns bool telling if self == other

        Returns True if the lines are the same type, and contain the same
        list of points.
        """
        if isinstance(other, self.__class__):
            # check to see if the end points are the same
            mypts = self.points
            otherpts = other.points
            same = True
            for (mypt, otherpt) in zip(mypts, otherpts):
                if mypt != otherpt:
                    same = False
                    break
            # flip to false if arc centers are different
            myactr = self.actr
            otheractr = other.actr
            if myactr != otheractr:
                same = False
            return same
        else:
            return False

    def __hash__(self):
        """Set the hash of the item to the id number.

        This allows one to make sets of these items.
        """
        return self.id

    def add_signline(self, signline):
        """Adds a signline to self.signlines."""
        if signline not in self.signlines:
            self.signlines.append(signline)

    def get_name(self):
        """Returns the arc name."""
        return 'L'+str(self.id)

    def length(self):
        """Returns the length of the arc."""
        a = self.pt(0)-self.actr
        b = self.pt(1)-self.actr
        ang = abs(a.ang_bet_rad(b))
        res = self.radius*ang
        return res

    def save_to_points(self):
        """This stores this line in the line's child points."""
        for point in self.allpoints:
            point.save_line(self)

    def set_ediv(self, ediv):
        """Sets the number of element divisions on the arc when meshing.

        Args:
            ediv (int): number of required elements on this arc
        """
        self.ediv = ediv

    def signed_copy(self, sign):
        """Returns a SignArc instance of this Arc with the passed sign."""
        return SignArc(self, sign)

    def pt(self, ind):
        """Returns the start or end point of the arc.

        Args:
            ind (int): index of the point we want, 0=start, 1=end

        Returns:
            Point: requested Point
        """
        return self.points[ind]

    def set_pt(self, ind, point):
        """Update the arc's point at ind to the new point.

        Args:
            ind (int): index of the point we're updating, 0=start, 1=end
            point (Point): new point we are assigning to Line
        """
        # remove this line from the old point
        self.points[ind].lines.remove(self)
        # set the new point
        self.points[ind] = point
        self.midpt = self.mid()
        self.save_to_points()
        other_point = self.pt(int(not ind))
        if point == other_point:
            other_point.lines.remove(self)
            for sline in self.signlines:
                sline.lineloop.remove(sline)
            if self.id != -1:
                if len(self.signlines) > 0:
                    area = self.signlines[0].lineloop.parent
                    area.part.fea.lines.remove(self)

    def mid(self):
        """Caclulates and returns the midpoint of the arc."""
        pt = self.get_pt_at(0.5)
        return pt

    def touches(self, other):
        """Checks if this arc is connected to a passed line.

        Args:
            other (Line or Arc): passed other line

        Returns:
            True: beginning or end point connect to other line
            False: lines are not connected
        """
        [p1, p2] = [self.pt(0), self.pt(1)]
        if p1 in other.points or p2 in other.points:
            return True
        else:
            return False

    @staticmethod
    def sgn(val):
        """Returns sign of the passed value.

        Helper function for self.interects code.

        Args:
            val (float): passed value

        Returns:
            1 or -1: sign of the passed value, assumes sign of zero = 1
        """
        if val < 0:
            return -1
        else:
            return 1

    def get_ang(self):
        """ Returns angle from beginning to end of arc from arc ctr in degrees

        Notes:
            Answer is between 180 and -180
            Positive is concave
            Negative is convex
        """
        a = self.pt(0)-self.actr
        b = self.pt(1)-self.actr
        ang = a.ang_bet_deg(b)
        return ang

    def get_ang_rad(self):
        """ Returns angle from beginning to end of arc from arc ctr in rad

        Notes:
            Answer is between pi and -pi
            Positive is concave
            Negative is convex
        """
        a = self.pt(0)-self.actr
        b = self.pt(1)-self.actr
        ang = a.ang_bet_rad(b)
        return ang

    def get_concavity(self, clockwise=True):
        """Returns concavity of this arc, 'concave' or 'convex'."""
        ang = self.get_ang()
        res = ''
        if clockwise:
            res = 'concave'
            if ang < 0:
                res = 'convex'
        else:
            res = 'convex'
            if ang < 0:
                res = 'concave'
        return res

    def get_pt_at(self, nondim):
        """Returns point on arc between start and end point at nondim dist.

        Args:
            nondim (float): distance between start and end point, 0=start, 1=end

        Returns:
            res (Point): requested point on arc
        """
        a = self.pt(0)-self.actr

        # neg is convex, pos is concave
        ang = self.get_ang()

        ang_delta = nondim*ang
        a.rot_ccw_deg(ang_delta)

        res = self.actr + a
        return res

    def get_perp_vec(self, pt):
        """ Returns vector perpendicular to current arc.

        Right hand rule is used:
            X+ is the direction of the arc from start to end.
            Z+ is out of the page
            Y+ is the direction of the perpendicular vector

        Args:
            pt (Point): point on arc where we want the perpendicular vector

        Returns:
            resv (Point): perpendicular vector
        """
        resv = None
        if self.concavity == 'convex':
            resv = pt - self.actr
        elif self.concavity == 'concave':
            resv = self.actr - pt
        return resv

    def get_tan_vec(self, point):
        """ Returns vector tangent to the current arc.

        Vector points from the passed point to a unit location further away.

        Args:
            point (Point): start point

        Returns:
            resv (Point): vector that is tangent to the current line
        """
        if self.concavity == 'concave':
            resv = self.actr - point
        elif self.concavity == 'convex':
            resv = point - self.actr
        if point == self.pt(0):
            resv.rot_ccw_deg(-90)
        else:
            resv.rot_ccw_deg(90)
        resv.make_unit()
        return resv

    def get_verts_codes(self, plot=True):
        """Returns a list of [verts, codes] vertices for plotting."""
        # http://www.spaceroots.org/documents/ellipse/node19.html
        # http://www.tinaja.com/glib/bezarc1.pdf
        # http://stackoverflow.com/questions/734076/geometrical-arc-to-bezier-curve

        # find the angle of rotation to the middle of the arc
        ang_offset = self.get_pt_at(0.5) - self.actr
        ang_offset = ang_offset.ang_deg()

        # store the half angle as theta
        theta = radians(abs(self.get_ang()))*0.5

        ax1 = (4 - cos(theta))/3
        rad1 = ((1 - cos(theta))*(cos(theta) - 3))/(3*sin(theta))
        pt1 = Point(rad1, ax1)*self.radius
        pt2 = Point(-rad1, ax1)*self.radius

        # fix the angles
        pt1.rot_ccw_deg(ang_offset)
        pt2.rot_ccw_deg(ang_offset)

        # fix the location
        pt1 += self.actr
        pt2 += self.actr

        # need to change the order if the arc is concave or convex
        verts = [pt1, pt2, self.pt(1)]
        if self.concavity == 'convex':
            verts = [pt2, pt1, self.pt(1)]

        if plot:
            verts = [point.axrad() for point in verts]
        else:
            verts = [point.radax() for point in verts]
        codes = [Path.CURVE4, Path.CURVE4, Path.CURVE4]
        return [verts, codes]

    def coincident(self, pt):
        """Checks to see if pt is on the arc.

        Args:
            pt (Point): input point to check

        Returns:
            bool: True if pt on arc or False if not
        """
        # (x - xc)^2 + (y - yc)^2 - r^2 = 0
        dist = pt - self.pt(0), pt - self.pt(1)
        dist = [item.length() for item in dist]
        if min(dist) < ACC:
            # True if point almost equal to end points
            return True
        remainder = abs((pt.x - self.actr.x)**2 + (pt.y - self.actr.y)**2
                        - self.radius**2)
        if remainder < ACC:
            # point is on the circle, but is it in the arc?
            # check if the angle traversed to the point is within the ang
            # traversed by the arc
            a1 = self.pt(0)-self.actr
            a2 = pt-self.actr
            ang_pt = a1.ang_bet_deg(a2)

            # this is the angle traversed byt the arc
            ang = self.get_ang()
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

    def intersects(self, other):
        """Checks if other line intersects this arc.

        Args:
            other (Line or Arc): other line to check

        Returns:
            Point or None:
                Point: intersection point
                None: None if lines don't intersect
        """
        if isinstance(other, Line) or isinstance(other, SignLine):
            #print('|C Checking arc-line intersection')
            #print(' %s' % self)
            # arc-line intersection
            # reset the math to be about the circle centroid
            # http://mathworld.wolfram.com/Circle-LineIntersection.html
            pt1 = other.pt(0) - self.actr
            pt2 = other.pt(1) - self.actr
            # for us x and y are reversed because of y axial
            (x1, y1, x2, y2) = (pt1.y, pt1.x, pt2.y, pt2.x)
            dx = x2 - x1
            dy = y2 - y1
            rad = self.radius
            dr = (dx**2 + dy**2)**(0.5)
            D = x1*y2 - x2*y1
            dist = [self.pt(0) - other.pt(0), self.pt(0) - other.pt(1),
                    self.pt(1) - other.pt(0), self.pt(1) - other.pt(1)]
            vals = [self.pt(0), self.pt(0), self.pt(1), self.pt(1)]
            dist = [item.length() for item in dist]
            for (distance, val) in zip(dist, vals):
                if distance < ACC:
                    return Point(val.x, val.y)
            discr = 1.0*(rad**2)*(dr**2) - D**2
            if discr < 0:
                # no intersection
                #print('none returned because of discriminant')
                return None
            else:
                # tangent or intersection
                rsq = rad**2
                drsq = dr**2
                x1 = (D*dy + self.sgn(dy)*dx*(rsq*drsq-D**2)**(0.5))/drsq
                y1 = (-D*dx + abs(dy)*(rsq*(dr**2)-D**2)**(0.5))/drsq
                x2 = (D*dy - self.sgn(dy)*dx*(rsq*drsq-D**2)**(0.5))/drsq
                y2 = (-D*dx - abs(dy)*(rsq*drsq-D**2)**(0.5))/drsq

                # convert result points back to global coordinates
                pt1 = Point(y1, x1) + self.actr
                pt2 = Point(y2, x2) + self.actr
                res = []
                # check that the resultant points are on line and arcs
                if self.coincident(pt1) and other.coincident(pt1):
                    for apoint in self.points:
                        dist = pt1 - apoint
                        dist = dist.length()
                        if dist < ACC:
                            pt1 = Point(apoint.x, apoint.y)
                            break
                    res.append(pt1)
                if pt1 != pt2 and self.coincident(pt2) and other.coincident(pt2):
                    for apoint in self.points:
                        dist = pt2 - apoint
                        dist = dist.length()
                        if dist < ACC:
                            pt2 = Point(apoint.x, apoint.y)
                            break
                    res.append(pt2)
                if len(res) == 1:
                    return res[0]
                elif len(res) > 1:
                    return res
                elif len(res) == 0:
                    #print('none returned because no point from math')
                    return None
        elif isinstance(other, Arc):
            # arc-arc intersection
            # I have not yet written the code here
            return None

    def __str__(self):
        """Returns string listing object type, name, and points"""
        name = self.get_name()
        start, end = self.pt(0), self.pt(1)
        actr = self.actr
        val = 'Arc %s, start: %s end: %s center: %s' % (name, start, end, actr)
        return val

class SignArc(Arc, base_classes.Idobj):
    """Makes a signed arc.

    Attributes:
        line (Arc): parent arc
        sign (int): 1 = positive, or -1 = negative
        lineloop (None or LineLoop): parent LineLoop, default is None
        nodes (list): nodes on the arc
        faces (list): list of element faces
        n1 (list): list of first order nodes on line
    """

    def __init__(self, parent_line, sign):
        self.line = parent_line
        self.sign = sign
        self.lineloop = None
        self.faces = []
        self.n1 = []
        base_classes.Idobj.__init__(self)

    def pt(self, ind):
        """Returns the point at the start or end of the line."""
        if self.sign == -1:
            return self.line.points[not ind]
        else:
            return self.line.points[ind]

    def set_pt(self, ind, point):
        """Sets the passed point to at the index location."""
        self.line.set_pt(ind, point)

    def set_ediv(self, ediv):
        """Apply the element divisions onto the parent line."""
        self.line.set_ediv(ediv)

    def set_lineloop(self, lineloop):
        """Sets the parent LineLoop"""
        # this is needed to cascade set ediv up to FEA model and down onto
        # other lines with this ediv
        self.lineloop = lineloop

    def get_name(self):
        """Returns the line name."""
        if self.sign == -1:
            return '-L'+str(self.line.id)
        else:
            return 'L'+str(self.line.id)

    def label(self, ax):
        """Labels the arc on the matplotlib ax axis.

        Args:
            ax (Matplotlib axis): Matplotlib axis
        """
        lvect = self.get_perp_vec(self.midpt)
        ang = lvect.ang_deg()
        horiz, vert = get_text_hv(ang)
        (axial, radial) = (self.midpt.y, self.midpt.x)
        ax.text(axial, radial, self.get_name(), ha=horiz, va=vert)

    def plot(self, ax, label=True):
        """Draws and labels the arc onto the passed Matplotlib axis.

        Args:
            ax: Matplotlib axis
            label (bool): if True, label the arc
        """
        ctr = self.actr
        vect1 = None
        sign = 1
        if self.concavity == 'concave':
            vect1 = self.pt(0)-ctr
        else:
            sign = -1
            vect1 = self.pt(1)-ctr
        rad = self.radius
        ang1 = (vect1.ang_deg())
        ang2 = ang1 + sign*self.get_ang()

        # matplotlib assumes ccw arc drawing, calculix assumes cw drawing
        a = AArc(xy=[ctr.y, ctr.x], width=2*rad, height=2*rad, angle=0,
                 theta1=ang1, theta2=ang2)
        ax.add_artist(a)
        if label:
            self.label(ax)

    def signed_copy(self, sign):
        """Returns a SignArc copy of this arc."""
        realsign = self.sign*sign
        return SignArc(self.line, realsign)

    @property
    def nodes(self):
        """Gets the line nodes."""
        return self.line.nodes

    @nodes.setter
    def nodes(self, nlist):
        """Sets the line nodes."""
        self.line.nodes = nlist

    @property
    def actr(self):
        """Return the arc center."""
        return self.line.actr

    @property
    def radius(self):
        """Return the arc radius."""
        return self.line.radius

    @property
    def midpt(self):
        """Return the mid point."""
        return self.line.midpt

    @property
    def edge(self):
        """Returns a bool saying if this is an edge line."""
        return self.line.edge

    @edge.setter
    def edge(self, isedge):
        """Sets the line nodes."""
        self.line.edge = isedge

    @property
    def points(self):
        """Return a list of the arc's points [pstart, pend]"""
        pstart = self.pt(0)
        pend = self.pt(1)
        return [pstart, pend]

    @property
    def concavity(self):
        """Return the arc's concavity 'concave' or 'convex'"""
        if self.sign == 1:
            return self.line.concavity
        else:
            if self.line.concavity == 'concave':
                return 'convex'
            else:
                return 'concave'

    def reverse(self):
        """Reverses self.line"""
        self.line.reverse()


class LineLoop(list):
    """Makes a LineLoop, which stores multiple Line/Arc or SignLine/SignArc.

    Args:
        items (list): Line or Arc or SignLine or SignArc
        hole_bool (bool): True if hole, False otherwise
        parent (Area or None): parent Area, defaults to None

    Attributes:
        closed (bool): True if closed
        area (float): the loop's area, can be pos or neg. Pos = cw, neg = ccw
        hole (bool): True if loop is a hole
        ccw (bool): True if loop is counter-clockwise
        center (Point): a point at the loop centroid
        parent (Area or None): the loop's parent Area
    """

    def __init__(self, items=[], hole_bool=False, parent=None):
        super().__init__(items)
        self.id = -1
        for item in self:
            if isinstance(item, SignLine) or isinstance(item, SignArc):
                item.lineloop = self
                if item not in item.line.signlines:
                    item.line.add_signline(item)
        self.hole = hole_bool
        self.parent = parent

    def set_id(self, id):
        """Sets the id number"""
        self.id = id

    def __hash__(self):
        """Sets the hash of the item to the id number.

        This allows one to make sets of these items.
        """
        return self.id

    @property
    def closed(self):
        """Returns True if Closed, False if not."""
        if len(self) < 2:
            return False
        else:
            for ind, curr_line in enumerate(self):
                lower_ind = ind-1
                prev_line = self[lower_ind]
                if curr_line.pt(0) != prev_line.pt(1):
                    return False
            return True

    @property
    def center(self):
        """Returns None or a Point at the line loop's center."""
        if self.closed == False:
            return None
        else:
            cx_sum = 0
            cy_sum = 0
            a_sum = 0
            arcs = []
            for item in self:
                if isinstance(item, Line) or isinstance(item, SignLine):
                    start, end = item.pt(0), item.pt(1)
                    det_val = det([[start.x, end.x], [start.y, end.y]])
                    cx_sum += (start.x + end.x)*det_val
                    cy_sum += (start.y + end.y)*det_val
                    a_sum += det_val
                elif isinstance(item, Arc) or isinstance(item, SignArc):
                    # concave: draw to actr, adder is neg
                    # convex:  draw to actr, adder is pos
                    tmp_list = [[item.pt(0), item.actr],
                                [item.actr, item.pt(1)]]
                    for [start, end] in tmp_list:
                        det_val = det([[start.x, end.x], [start.y, end.y]])
                        cx_sum += (start.x + end.x)*det_val
                        cy_sum += (start.y + end.y)*det_val
                        a_sum += det_val
                    arcs.append(item)
            area_val = a_sum*0.5
            val_list = []
            if area_val != 0:
                cx_val = (1/(6*area_val))*cx_sum
                cy_val = (1/(6*area_val))*cy_sum
                val_list.append([area_val, cx_val, cy_val])
            for arc in arcs:
                ang = arc.get_ang_rad()
                area_val = -0.5*ang*(arc.radius**2)
                half_ang = abs(ang)*0.5
                vect = arc.midpt - arc.actr
                vect.make_unit()
                mag = 2*arc.radius*sin(half_ang)/(3*half_ang)
                vect = vect*mag
                center = arc.actr + vect
                val_list.append([area_val, center.x, center.y])
            a_sum = sum([aval for [aval, cx, cy] in val_list])
            cxa_sum = sum([cx*aval for [aval, cx, cy] in val_list])
            cya_sum = sum([cy*aval for [aval, cx, cy] in val_list])
            cx_val = cxa_sum/a_sum
            cy_val = cya_sum/a_sum
            return Point(cx_val, cy_val)

    @property
    def ccw(self):
        """Returns bool telling if the area is closed and clockwise."""
        if self.closed == False:
            return False
        else:
            area_val = self.area
            if area_val < 0:
                return True
            else:
                return False

    @property
    def area(self):
        """Returns the loop's signed area.

        Returns 0 if not closed.
        Pos = cw, neg = ccw

        Formula from: http://mathworld.wolfram.com/PolygonArea.html
        """
        res = 0
        adder = 0
        arcs = []
        if self.closed == False:
            return res
        else:
            for item in self:
                if isinstance(item, Line) or isinstance(item, SignLine):
                    start, end = item.pt(0), item.pt(1)
                    det_matrix = [[start.x, end.x], [start.y, end.y]]
                    res += det(det_matrix)
                elif isinstance(item, Arc) or isinstance(item, SignArc):
                    # concave: draw to actr, adder is neg
                    # convex:  draw to actr, adder is pos
                    start, end = item.pt(0), item.actr
                    det_matrix = [[start.x, end.x], [start.y, end.y]]
                    res += det(det_matrix)
                    start, end = item.actr, item.pt(1)
                    det_matrix = [[start.x, end.x], [start.y, end.y]]
                    res += det(det_matrix)
                    arcs.append(item)
            res = 0.5*res
            # need to loop through arcs here, now that we know if cw or ccw
            for arc in arcs:
                ang = arc.get_ang_rad()
                tmp_adder = -0.5*ang*(arc.radius**2)
                adder += tmp_adder
            res = res + adder
            return res

    def reverse(self):
        """Reverse the loop direction cw<-->ccw

        This modifies the order of lines, and flips each line.
        """
        for line in self:
            line.reverse()
        super().reverse()

    def set_parent(self, parent):
        """Sets the line loop's parent part."""
        self.parent = parent

    def get_patch(self):
        """Returns a patch of the LineLoop, or None of not closed."""
        if self.closed:
            codes, verts = [], []
            # start poly
            codes += [Path.MOVETO]
            verts += [self[0].pt(0).radax()]
            for line in self:
                if isinstance(line, SignLine) or isinstance(line, Line):
                    verts.append(line.pt(1).radax())
                    codes.append(Path.LINETO)
                elif isinstance(line, SignArc) or isinstance(line, Arc):
                    lverts, lcodes = line.get_verts_codes(plot=False)
                    verts += lverts
                    codes += lcodes
                if line == self[-1]:
                    # close poly if at last line
                    verts.append((0, 0))
                    codes.append(Path.CLOSEPOLY)
            path = Path(verts, codes)
            patch = PathPatch(path, facecolor='yellow', linewidth=0.0)
            return patch
        return None

    def contains_point(self, point):
        """Returns bool telling if pt is inside the LineLoop.

        Args:
            pt (Point): point we are checking

        Returns:
            contains (bool): True if in this area False if not
        """
        contains = False
        if self.closed:
            patch = self.get_patch()
            contains = patch.contains_point(point.radax())
        return contains

    def inside(self, other):
        """Checks to see if self is inside other.

        Only checks self's points.
        """
        if self.closed == False:
            return False
        if isinstance(other, LineLoop):
            points = set()
            for sline in self:
                points.update(sline.points)
            patch = other.get_patch()
            for point in points:
                contains = patch.contains_point(point.radax())
                if contains == False:
                    return False
            return True
        else:
            return False

    def append(self, item):
        """Adds an item to this LineLoop

        Args:
          item (SingArc or SignLine): item to add to the list
        """
        if isinstance(item, SignLine) or isinstance(item, SignArc):
            item.set_lineloop(self)
        super().append(item)
        return item

    def remove(self, item):
        """Removes an SignLine or SignArc item from this LineLoop

        Args:
          item (SignArc or SignLine): item to remove
        """
        if isinstance(item, SignLine) or isinstance(item, SignArc):
            item.set_lineloop(None)
            if item not in item.line.signlines:
                print('''===THIS SLINE NOT IN SLINE.LINE.SIGNLINES!===''')
                print('Removing sline %s from loop in area %s'
                      % (item.get_name(), self.parent.get_name()))
                print('SLINE.LINE %s' % item.line)
                print('SLINE.LINE.SIGNLINES %r' % item.line.signlines)
        super().remove(item)

    def insert(self, index, item):
        """Inserts a SignLine or SignArc item into this LineLoop

        Args:
          item (SingArc or SignLine): item to insert
        """
        if isinstance(item, SignLine) or isinstance(item, SignArc):
            item.set_lineloop(self)
        super().insert(index, item)

    def __str__(self):
        """Prints the LineLoop definition."""
        res = 'Number of items %i\n' % len(self)
        for line in self:
            res += '%s\n' % line
        return res

class Area(base_classes.Idobj):
    """ Makes an area.

    Area is closed in a clockwise direction.

    Args:
        part (Part): parent Part
        line_list (list): list of Lines and Arcs that close the area

    Attributes:
        part (Part): parent part
        closed (boolean): True if closed, False if not
        exlines (LineLoop): LineLoop of signlines that define the area exterior
        signlines (list): list of signed lines or arcs that define the area
        lines (list): a list of all the lines that make the area, includes hole
            lines
        point (list): a list of all points making the area, excludes arc centers
        allpoints (list): a list of all points, includes arc centers
        holepoints (list): a list of hole points, excludes arc centers
        matl (Matl): the material fo the area
        etype (str): element type of area. Options are:
            'plstress': plane stress
            'plstrain': plane strain
            'axisym': axisymmetric
        area (double): in-plane part area (ballpark number, arcs not calc right)
        holes (list): list of LineLoop of signlines
        center (Point or None): the area centroid point. None before closed.
        nodes (list): child mesh nodes that are in the area
        elements (list): child elements that are in the area
    """

    def __init__(self, part, line_list=[]):
        base_classes.Idobj.__init__(self)
        self.part = part # parent part reference
        self.closed = False
        self.matl = None
        self.etype = None
        self.area = 0.0
        self.holes = []
        self.center = None
        self.nodes = []
        self.elements = []
        exloop = LineLoop(line_list, False, self)
        self.part.fea.register(exloop)
        self.exlines = exloop
        if self.exlines.closed == True:
            self.close()

    def __hash__(self):
        """Returns the item's id as its hash."""
        return self.id

    @property
    def lines(self):
        """Returns a list of all Line and Arc under this area."""
        mylist = []
        loops = [self.exlines] + self.holes
        for loop in loops:
            line_list = [line.line for line in loop]
            mylist += line_list
        return mylist

    @property
    def signlines(self):
        """Returns a list of all SignLine and SignArc under this area."""
        mylist = []
        loops = [self.exlines] + self.holes
        for loop in loops:
            mylist += loop
        return mylist

    @property
    def points(self):
        """Returns a list of all points under this area; excludes arc centers."""
        pts = set()
        for line in self.lines:
            pts.update(line.points)
        return list(pts)

    @property
    def allpoints(self):
        """Returns a list of all points under this area; includes arc centers."""
        pts = set()
        for line in self.lines:
            pts.update(line.allpoints)
        return list(pts)

    @property
    def holepoints(self):
        """Returns a list of all points in this area's holes."""
        pts = set()
        for hole in self.holes:
            for sline in hole:
                pts.update(sline.points)
        return list(pts)

    def __str__(self):
        """Returns a string listing the area's name."""
        return self.get_name()

    def get_name(self):
        """Returns the area name."""
        return 'A'+str(self.id)

    def close(self):
        """Closes the area."""
        # sets closed=True
        # calculates the area center
        self.closed = True
        if self.exlines.ccw == True:
            self.exlines.reverse()
        self.area, self.center = self.calc_area_center()

    def line_from_startpt(self, pt):
        """Returns the signline in the area that start at the passed pt.

        Args:
            pt (Point): the start point of the line one wants

        Returns:
            match (None or SignLine or SignArc):
            |  Returns None if no line is found, otherwise the correct line or
            |  arc which starts with pt will be
            |  returned.
        """
        # returns the line that starts on the given point
        matches = []
        area_slines = self.signlines
        for sline in area_slines:
            if sline.pt(0) == pt:
                matches.append(sline)
        if len(matches) == 1:
            return matches[0]
        elif len(matches) > 1:
            # we have a saw cut in the area and we're cutting on one of its pts
            print('>1 slines starting with %s' % pt.get_name())
            for match in matches:
                slines = set(match.line.signlines)
                if len(slines) == 2 and slines.issubset(area_slines):
                    # this is the saw cut one
                    pass
                else:
                    # we want the non saw cut one
                    return match
        return None

    def calc_area_center(self):
        """Calculates and returns the area and centroid Point.

        Returns:
            list: [area, Point]
        """
        loops = [self.exlines] + self.holes
        val_list = []
        for loop in loops:
            if loop.closed == True:
                aval, cval = loop.area, loop.center
                val_list.append([aval, cval])
        a_sum = sum([aval[0] for aval in val_list])
        cxa_sum = sum([center.x*aval for [aval, center] in val_list])
        cya_sum = sum([center.y*aval for [aval, center] in val_list])
        cx_val = cxa_sum/a_sum
        cy_val = cya_sum/a_sum
        center = Point(cx_val, cy_val)
        return [a_sum, center]

    def get_maxlength(self):
        """Returns max distance between points in an area."""
        pts = self.points
        maxlen = 0.0
        # loop through points checking dist to next point
        for ind, p1 in enumerate(pts[:-1]):
            for p2 in pts[ind:]:
                vect = p1 - p2
                dist = vect.length()
                if dist > maxlen:
                    maxlen = dist
        return maxlen

    def rem_sline(self, sline):
        """Removes a sline from this area."""
        loops = [self.exlines] + self.holes
        for loop in loops:
            if sline in loop:
                print('Line %s removed from %s'
                      % (sline.get_name(), self.get_name()))
                loop.remove(sline)
        if sline.id != -1:
            line = sline.line
            self.part.fea.signlines.remove(sline)
            if len(line.signlines) == 0:
                self.part.fea.lines.remove(line)

    def add_sline(self, signline):
        """Adds signline to the area definition.

        SignLine is added to the end of the lines list.
        Area is closed if needed.
        """
        # this adds a line to the area
        self.exlines.append(signline)
        if self.exlines.closed:
            self.close()

    def add_hole_sline(self, signline):
        """Adds signline to the area hole definition.

        Line is added to the end of the lines list.
        Area is closed if needed.

        Returns:
            closed (bool): boolean telling if the hole was closed
        """
        make_new = (len(self.holes) == 0)
        if make_new == False:
            if self.holes[-1].closed == True:
                make_new = True
        if make_new:
            hole = LineLoop([], hole_bool=True, parent=self)
            hole = self.part.fea.lineloops.append(hole)
            self.holes.append(hole)
        self.holes[-1].append(signline)
        # reverse the hole if it's closed in the cw direction
        if self.holes[-1].closed == True:
            if self.holes[-1].ccw == False:
                self.holes[-1].reverse()
            self.area, self.center = self.calc_area_center()
        return self.holes[-1].closed

    def get_patch(self):
        """Returns a patch of the area."""
        if self.closed:
            loops = [self.exlines] + self.holes
            codes, verts = [], []
            for loop in loops:
                # start poly
                codes += [Path.MOVETO]
                verts += [loop[0].pt(0).axrad()]
                for sline in loop:
                    if isinstance(sline, SignLine):
                        verts.append(sline.pt(1).axrad())
                        codes.append(Path.LINETO)
                    elif isinstance(sline, SignArc):
                        lverts, lcodes = sline.get_verts_codes()
                        verts += lverts
                        codes += lcodes
                    if sline == loop[-1]:
                        # close poly if at last line
                        verts.append((0, 0))
                        codes.append(Path.CLOSEPOLY)
            path = Path(verts, codes, _interpolation_steps=2)
            patch = PathPatch(path, linewidth=0.0)
            return patch

    def label(self, axis):
        """Labels the area on a Matplotlib axis

        Args:
            axis (Matplotlib Axis): Matplotlib Axis
        """
        axis.text(self.center.y, self.center.x, self.get_name(),
                  ha='center', va='center')

    def plot(self, axis, label=True, color='yellow'):
        """Plots the area on a Matplotlib axis.

        Args:
            axis (Matplotlib axis): Matplotlib axis
            label (bool): if True, label area
        """
        # http://stackoverflow.com/questions/8919719/how-to-plot-a-complex-polygon
        if self.closed:
            # plot area patch
            patch = self.get_patch()
            patch.set_color(color)
            axis.add_patch(patch)
            # label ares
            if label:
                self.label(axis)

    def contains_point(self, point):
        """Returns bool telling if pt is inside the area.

        Args:
            pt (Point): point we are checking

        Returns:
            contains (bool): True if in this area False if not
        """
        patch = self.get_patch()
        contains = patch.contains_point(point.axrad())
        return contains

    def line_insert(self, lgiven, lnew, after=True):
        """Inserts line lnew before or after the given line.

        Args:
            lgiven (SignLine or SignArc): the line we are inserting relative to
            lnew (SignLine or SignArc): the new line we will insert
            after (bool): if True insert lnew after lgiven, false insert before

        Returns:
            bool: True if succeeded, False if failed
        """
        if lgiven in self.signlines:
            print('Inserting line %s into area %s next to %s' %
                  (lnew.get_name(), self.get_name(), lgiven.get_name()))
            lineloop = lgiven.lineloop
            type = 'exterior lines'
            if lineloop.hole:
                type = 'hole lines'
            ind = lineloop.index(lgiven)
            if after:
                ind += 1
            print(' Inserting line into %s outer lines' % self.get_name())
            lineloop.insert(ind, lnew)
            return True
        else:
            print('ERROR: Passed line was not in %s!' % self.get_name())
            print('You must pass a line in this area!')
            return False

    def set_child_ccxtypes(self):
        """Assigns the child element ccx types."""
        if self.etype != None and self.elements != []:
            # if etype is set and elements have been assigned
            for element in self.elements:
                element.set_ccxtype(self.etype)

    def set_etype(self, etype):
        """Sets the area etype.

        Sets the ccxtype on elements if they have been assigned to this area.

        Args:
            etype (str): element type
                'plstress' plane stress
                'plstrain': plane strain
                'axisym': axisymmetric
        """
        self.etype = etype
        self.set_child_ccxtypes()

    def update(self, line_list):
        """Updates the area's external lines to the passed list.

        If the area is cosed, self.close() will be called.

        Args:
            line_list (list): list of SignLine and SignArc items
        """
        to_rem = list(self.exlines)
        for sline in to_rem:
            self.exlines.remove(sline)
        for sline in line_list:
            self.exlines.append(sline)
        if self.exlines.closed == True:
            self.close()
