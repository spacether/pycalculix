"""This module stores geometry classes, which are used to make parts.
"""
from math import atan2, pi, cos, sin, radians
from matplotlib.patches import Arc as AArc
from matplotlib.patches import PathPatch
from matplotlib.path import Path

from . import base_classes

#accuracy for small numbers in math below
ACC = .00001

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
    """

    def __init__(self, x, y, z=0):
        self.x = x
        self.y = y
        self.z = z
        self.nodes = []
        self.lines = []
        self.arc_center = False
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
        """Saves the line in th epoint."""
        if line not in self.lines:
            self.lines.append(line)

    def length(self):
        """Returns the length of this point or vector."""
        res = (self.x**2 + self.y**2)**0.5
        return res

    def label(self, ax):
        """Labels the point on a Matplotlib axis.

        Args:
            ax (Matplotlib Axis): Matplotlib Axis
        """
        ax.annotate(self.get_name(), (self.y, self.x))

    def plot(self, ax, label=True):
        """Plots the point on the passed matplotlib axis.

        Args:
            ax (Matplotlib axis): plate to plot the point
            pnum (bool): True turns on point labeling
        """
        ax.scatter(self.y, self.x)
        if label:
            self.label(ax)

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

    def rot_ccw_deg(self, ang):
        """Rotates the current vector by ccw degrees about 0,0. Returns self.

        Args:
            ang (float): angle or rotation in degrees, counter-clockwise (ccw)
        """
        ang = radians(ang)
        ax = self.y*cos(ang) - self.x*sin(ang)
        rad = self.y*sin(ang) + self.x*cos(ang)
        ax = round(ax, 5)
        rad = round(rad, 5)
        (self.x, self.y) = (rad, ax)
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
        val = 'Point, id %i, (x,y)=(%f,%f)' % (self.id, self.x, self.y)
        return val


class Line(base_classes.Idobj):
    """Stores a line from points p1 to p2.

    Args:
        p1 (Point): first point
        p2 (Point): second point
        sign (int): 1 or -1, defines direction of line

            - +1: p1 -> p2
            - -1: p2 -> p1

    Attributes:
        points (list of Point) : list of the line's points [p1,p2]
        allpoints (list of point): same as points
        midpt (Point): a mid-point between p1 and p2
        nodes (list of Node): a list of meshed nodes on the line
        faces (list of Face): a list of meshed faces on the line
    """

    def __init__(self, p1, p2):
        self.points = [p1, p2]
        self.midpt = self.mid()
        self.ediv = None
        self.nodes = []
        self.edge = True
        self.signlines = []
        base_classes.Idobj.__init__(self)

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

    def set_pt(self, ind, pt):
        """Update the line's point at ind to the new pt.

        Args:
            ind (int): index of the point we're updating, 0=start, 1=end
            pt (Point): new point we are assigning to Line
        """
        # used to update specific points in the line
        self.points[ind] = pt
        self.midpt = self.mid()

    def mid(self):
        """Caclulates and returns a midpoint between start and end points."""
        pt = (self.pt(0) + self.pt(1))*0.5
        return pt

    def touches(self, other):
        """Checks if a line is connected to a passed line.

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
        p0 = self.pt(0) + lvect
        p1 = self.pt(1) + lvect
        tmpline = Line(p0, p1)
        return tmpline

    def get_abc(self):
        """ Returns a list of abc terms for a line definition.

        ax + by + c = 0

        Returns:
            [a, b, c]: list of the above line terms
        """
        lpt = self.pt(1)-self.pt(0)
        dx = lpt.x
        dy = lpt.y
        (a, b, c) = (0.0, 0.0, 0.0)
        if dy == 0:
            # vertical radial line
            (b, c) = (1.0, self.pt(0).y*-1.0)
        elif dx == 0:
            # horizontal axial line
            (a, c) = (1.0, self.pt(0).x*-1.0)
        else:
            slope = dy/dx
            offset = self.pt(0).y - slope*self.pt(0).x
            (a, b, c) = (1.0, -1.0*slope, -1.0*offset)
        return [a, b, c]

    def get_perp_vec(self, pt=None):
        """ Returns vector perpendicular to current line.

        Vector is created by rotating the current line ccw 90 degrees
        around the start point.

        Args:
            pt (Point): defaults to None, unused for Line, used for Arc

        Returns:
            lvect (Point): vector that is perpendiculat to the current line
        """
        lvect = self.pt(1) - self.pt(0)
        lvect.rot_ccw_deg(90)
        return lvect

    def arc_tang_intersection(self, pt, mag):
        """Returns an intersection point on this line of an arc centered at pt

        Args:
            pt (Point): arc center point
            mag (float): passed radius distance from the arc center to the line

        Returns:
            newpt (Point or None):
                Point: the intersection point
                None: returns None if no intersection point exists
        """
        v = self.get_perp_vec()
        v.make_unit()
        p2 = pt+v*(-2*mag)

        tmpline = Line(pt, p2)
        newpt = tmpline.intersects(self)
        if type(newpt) == type(None):
            print('Intersection failed!')
        return newpt

    def coincident(self, pt):
        """Checks to see if pt is on the line.

        Args:
            pt (Point): input point to check

        Returns:
            bool: True if pt on line or False if not
        """
        [a, b, c] = self.get_abc()
        remainder = abs(a*pt.x + b*pt.y + c)
        if remainder < ACC:
            # point is on the line equation, but is it between end pts
            lvect = self.pt(1) - self.pt(0)
            # point0 + lvect*term = pt
            nondim = pt - self.pt(0)
            nondim = nondim/lvect
            if 0 <= nondim <= 1.0:
                # we're on the line
                return True
            else:
                # we're off of the line
                return False
        else:
            # point is not on line equation
            return False

    def intersects(self, other):
        """Checks if other line intersects this line.

        Args:
            other (Line or Arc): other line to check

        Returns:
            Point or None:
                Point: intersection point
                None: None if lines don't intersect
        """
        # returns intersection point between this line and other line
        # other is a straight line
        if isinstance(other, Line) or isinstance(other, SignLine):
            slope1 = self.get_abc()
            slope2 = other.get_abc()
            if slope1[:2] == slope2[:2]:
                # lines are paralel and will not have an intersection point
                return None
            else:
                # parametric intersection calculation
                # p1 -> p2 = ta parametric, p2 -> p3 tb parametric
                (p1, p2, p3, p4) = (self.pt(0), self.pt(1), other.pt(0), other.pt(1))
                (x1, x2, x3, x4) = (p1.y, p2.y, p3.y, p4.y)
                (y1, y2, y3, y4) = (p1.x, p2.x, p3.x, p4.x)
                (x31, x21, x43) = (x3-x1, x2-x1, x4-x3)
                (y31, y21, y43) = (y3-y1, y2-y1, y4-y3)
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

    def __str__(self):
        """Returns string listing object type, id number, and points"""
        p0 = self.pt(0)
        p1 = self.pt(1)
        val = 'Line %s  p0: %s p1: %s' % (self.get_name(), p0, p1)
        return val

class SignLine(Line, base_classes.Idobj):
    """Makes a signed line.

    Attributes:
        line (Line): parent line
        sign (int): 1 = positive, or -1 = negative
        parent (None or Area): parent Area, start out as None
        faces (list): list of element faces
        n1 (list): list of first order nodes on line
        nodes (list): nodes on the line
    """
    def __init__(self, parent_line, sign):
        self.line = parent_line
        self.sign = sign
        self.parent = None
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

    def set_parent(self, p):
        """Sets area parent."""
        # this is needed to cascade set ediv up to FEA model and down onto
        # other lines with this ediv
        self.parent = p

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

class Arc(base_classes.Idobj):
    """Makes an arc from points p1 to p2 about center actr.

    Arcs should be 90 degrees maximum.

    Arguments:
        p1 (Point): first point
        p2 (Point): second point
        actr (Point): arc center
        sign (int): 1 or -1, defines direction of line

            - +1: p1 -> p2
            - -1: p2 -> p1

    Attributes:
        points (list of Point) : list of the line's points [p1,p2]
        allpoints (list of Point): [p1, p2, actr]
        sign (int): 1 or -1, defines direction of line

            - +1: p1 -> p2
            - -1: p2 -> p1

        actr (Point): arc center
        radius (float): radius of the arc
        concavity (str):

            - 'concave': concave arc
            - 'convex': convex arc

        midpt (Point): a mid-point between p1 and p2
        nodes (list of Node): a list of meshed nodes on the line
        faces (list of Face): a list of meshed faces on the line
    """

    def __init__(self, p1, p2, actr):
        self.points = [p1, p2]
        self.actr = actr
        self.actr.arc_center = True
        self.radius = (p1-actr).length()
        self.concavity = self.get_concavity()
        self.midpt = self.mid()
        self.ediv = None
        self.nodes = []
        self.edge = True
        self.signlines = []
        base_classes.Idobj.__init__(self)

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
        ang = a.ang_bet_rad(b)
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

    def get_concavity(self):
        """Returns concavity of this arc, 'concave' or 'convex'."""
        ang = self.get_ang()
        res = 'concave'
        if ang < 0:
            res = 'convex'
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

    def coincident(self, pt):
        """Checks to see if pt is on the arc.

        Args:
            pt (Point): input point to check

        Returns:
            bool: True if pt on arc or False if not
        """
        # (x - xc)^2 + (y - yc)^2 - r^2 = 0
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
        if isinstance(other, Line):
            # arc-line intersection
            # reset the math to be about the circle centroid
            # formula from: http://mathworld.wolfram.com/Circle-LineIntersection.html
            p1 = other.pt(0) - self.actr
            p2 = other.pt(1) - self.actr
            # for us x and y are reversed because of y axial
            (x1, y1, x2, y2) = (p1.y, p1.x, p2.y, p2.x)
            dx = x2 - x1
            dy = y2 - y1
            r = self.radius
            dr = (dx**2 + dy**2)**(0.5)
            D = x1*y2 - x2*y1
            discr = 1.0*(r**2)*(dr**2) - D**2
            if discr < 0:
                # no intersection
                return None
            else:
                # tangent or intersection
                x1 = (D*dy+self.sgn(dy)*dx*((r**2)*(dr**2)-D**2)**(0.5))/(dr**2)
                y1 = (-D*dx+abs(dy)*dx*((r**2)*(dr**2)-D**2)**(0.5))/(dr**2)
                x2 = (D*dy-self.sgn(dy)*dx*((r**2)*(dr**2)-D**2)**(0.5))/(dr**2)
                y2 = (-D*dx-abs(dy)*dx*((r**2)*(dr**2)-D**2)**(0.5))/(dr**2)

                # convert result points back to global coordinates
                p1 = Point(y1, x1) + self.actr
                p2 = Point(y2, x2) + self.actr
                res = []
                # check that the resultant points are on line and arcs
                if self.coincident(p1) and other.coincident(p1):
                    res.append(p1)
                if p1 != p2 and self.coincident(p2) and other.coincident(p2):
                    res.append(p2)
                if len(res) == 1:
                    return res[0]
                elif len(res) > 1:
                    return res
                elif len(res) == 0:
                    return None
        elif isinstance(other, Arc):
            # arc-arc intersection
            # I have not yet written the code here
            return None

    def __str__(self):
        """Returns string listing object type, id number, and points"""
        p0 = self.pt(0)
        p1 = self.pt(1)
        actr = self.actr
        val = 'Arc, id %i p0: %s p1: %s actr: %s' % (self.id, p0, p1, actr)
        return val

class SignArc(Arc, base_classes.Idobj):
    """This class is for a signed line."""

    def __init__(self, parent_line, sign):
        self.line = parent_line
        self.sign = sign
        self.parent = None
        self.nodes = []
        self.faces = []
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

    def set_parent(self, p):
        """Sets area parent."""
        # this is needed to cascade set ediv up to FEA model and down onto
        # other lines with this ediv
        self.parent = p

    def get_name(self):
        """Returns the line name."""
        if self.sign == -1:
            return '-L'+str(self.line.id)
        else:
            return 'L'+str(self.line.id)

    def get_verts_codes(self):
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
        p1 = Point(rad1, ax1)*self.radius
        p2 = Point(-rad1, ax1)*self.radius

        # fix the angles
        p1.rot_ccw_deg(ang_offset)
        p2.rot_ccw_deg(ang_offset)

        # fix the location
        p1 += self.actr
        p2 += self.actr
        verts = [p1.axrad(), p2.axrad(), self.pt(1).axrad()]
        codes = [Path.CURVE4, Path.CURVE4, Path.CURVE4]
        return [verts, codes]

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


class Lineloop(list):
    """Makes a line loop, which stores multiple Line/Arc or SignLine/SignArc.
    """
    
    def __init__(self, items, hole_bool=False):
        list.__init__(self, items)
        self.closed = False
        self.area = 0
        self.hole = hole_bool
    
    @property
    def closed(self):
        """Returns True if Closed, False if not."""
        print('ran closure function')
        if len(self) < 2:
            return False
        else:
            for ind, curr_line in enumerate(self):
                lower_ind = ind-1
                prev_line = self[lower_ind]
                if curr_line.pt(0) != prev_line.pt(1):
                    return False
            return True

    def signed_area(self):
        """Returns the loop's signed area.
        
        Returns 0 if not closed.
        """
        if self.closed == False:
            return 0
        else:
            pass

class Area(base_classes.Idobj):
    """ Makes an area.

    Area is closed in a clockwise direction.

    Args:
        p (Part): parent part
        line_list (list): list of Lines and Arcs that close the area

    Attributes:
        p (Part): parent part
        closed (boolean): True if closed, False if not
        exlines (list): list of signlines that define the area exterior
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
        holes (list): list of list of signlines
        center (Point or None): the area centroid point. None before closed.
        nodes (list): child mesh nodes that are in the area
        elements (list): child elements that are in the area
    """

    def __init__(self, parent, line_list=[]):
        # initialtes the area
        base_classes.Idobj.__init__(self)
        self.parent = parent # parent part reference
        self.closed = False
        self.exlines = []
        self.matl = None
        self.etype = None
        self.area = 0.0
        self.holes = []
        self.center = None
        self.nodes = []
        self.elements = []
        self.update(line_list)

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
        self.calc_center()

    def line_from_startpt(self, pt):
        """Returns the signline in the area that start at the passed pt.

        Args:
            pt (Point): the start point of the line one wants

        Returns:
            None or line (SignLine or SignArc): Returns None if no line is found
                , otherwise the correct line or arc which starts with pt will be
                returned.
        """
        # returns the line that starts on the given point
        for line in self.signlines:
            if line.pt(0) == pt:
                return line
        return None

    def calc_center(self):
        """Sets self.center the area centroid Point.

        Attributes area and centroid are set by this method

        area (double): in-plane part area
        center (Point): the area centroid point
        """

        # loop to get area and centroid
        cx = 0.0
        cy = 0.0
        area = 0.0
        for line in self.exlines:
            pt1 = line.pt(0)
            pt2 = line.pt(1)
            term1 = (pt1.x*pt2.y)-(pt2.x*pt1.y)
            area += 1.0*term1
            cx += (pt1.x + pt2.x)*term1
            cy += (pt1.y + pt2.y)*term1
        area = 0.5*area
        cx = (1/(6.0*area))*cx
        cy = (1/(6.0*area))*cy
        self.area = area
        self.center = Point(cx, cy)

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

    def add_sline(self, signline):
        """Adds signline to the area definition.

        SignLine is added to the end of the lines list.
        Area is closed if needed.
        """
        # this adds a line to the area
        signline.set_parent(self)
        self.exlines.append(signline)
        # check to see if the area is closed
        if signline.pt(1) == self.exlines[0].pt(0):
            self.close()

    def add_hole_sline(self, signline):
        """Adds signline to the area hole definition.

        Line is added to the end of the lines list.
        Area is closed if needed.

        Returns:
            closed (bool): boolean telling if the hole was closed
        """
        signline.set_parent(self)
        if len(self.holes) == 0:
            self.holes.append([])
        self.holes[-1].append(signline)
        closed = False
        print('Hole line end point %s' % signline.pt(1))
        print('Hole first point: %s' % self.holes[-1][0].pt(0))
        if signline.pt(1) == self.holes[-1][0].pt(0):
            # last point == first point
            closed = True
        print('Equality: %s' % str(closed))
        return closed

    def get_patch(self):
        """Returns a patch of the area."""
        if self.closed:
            loops = [self.exlines] + self.holes
            codes, verts = [], []
            for loop in loops:
                # start poly
                codes += [Path.MOVETO]
                verts += [loop[0].pt(0).axrad()]
                for l in loop:
                    if isinstance(l, SignLine):
                        verts.append(l.pt(1).axrad())
                        codes.append(Path.LINETO)
                    elif isinstance(l, SignArc):
                        lverts, lcodes = l.get_verts_codes()
                        verts += lverts
                        codes += lcodes
                    if l == loop[-1]:
                        # close poly if at last line
                        verts.append((0, 0))
                        codes.append(Path.CLOSEPOLY)
            path = Path(verts, codes)
            patch = PathPatch(path, facecolor='yellow', linewidth=0.0)
            return patch

    def plot(self, ax, label=True, fill=True):
        """Plots the area on a Matplotlib axis.

        Args:
            ax (Matplotlib axis): Matplotlib axis
            label (bool): if True, label area
            fill (bool): if True, fill in area
        """
        # http://stackoverflow.com/questions/8919719/how-to-plot-a-complex-polygon
        if self.closed:
            # plot area patch
            if fill:
                patch = self.get_patch()
                ax.add_patch(patch)
            # label ares
            if label:
                ax.text(self.center.y, self.center.x, self.get_name(),
                        ha='center', va='center')

    def contains_point(self, point):
        """Returns bool telling if pt is inside the area.

        Args:
            pt (Point): point we are checking

        Returns:
            contains (bool): True if in this area False if not
        """
        patch = self.get_patch()
        ax = point.y
        rad = point.x
        contains = patch.contains_point((ax, rad))
        return contains

    def line_insert(self, lgiven, lnew, after=True):
        """Inserts line lnew before or after the given line.

        Args:
            lgiven (SignLine or SignArc): the line we are inserting relative to
            lnew (Line or Arc): the new line we will insert
            after (bool): if True insert lnew after lgiven, false insert before

        Returns:
            bool: True if succeeded, False if failed
        """
        # this inserts a line after the given line
        if lgiven in self.exlines:
            lnew.set_parent(self)
            ind = self.exlines.index(lgiven)
            if after:
                ind += 1
            self.exlines.insert(ind, lnew)
            return True
        else:
            for hole in self.holes:
                if lgiven in hole:
                    lnew.set_parent(self)
                    ind = hole.index(lgiven)
                    if after:
                        ind += 1
                    hole.insert(ind, lnew)
                    return True
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
        """Updates the area's signlines list to the passed list.

        The area's signlines list is replaced with the passed list.
        If the area is cosed, self.close() will be called.

        Args:
            line_list (list): list of SignLine and SignArc items
        """
        # adds the line list and checks for closure
        self.exlines = line_list
        for sline in self.exlines:
            sline.set_parent(self)
        if len(line_list) > 2:
            if self.exlines[0].pt(0) == self.exlines[-1].pt(1):
                self.close()
