from math import atan2, pi, cos, sin, radians
from matplotlib.patches import Arc as AArc

from . import base_classes

#accuracy for small numbers in math below
_acc = .00001

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
        list length is 1.
    """
    
    def __init__(self, x, y, z=0):
        self.x = x
        self.y = y
        self.z = z
        self.nodes= []
        base_classes.Idobj.__init__(self)
        
    def get_name(self):
        """Returns item name."""
        return 'P'+str(self.id)
    
    def __eq__(self, other):
        """Checks x and y equality to other point.
        
        Args:
          other: passed point        
        """        
        if self.x == other.x and self.y == other.y:
            return True
        return False
        
    def __add__(self, other):
        """Returns new point = self + other.
        
        Args:
          other: passed point
        
        Returns:
          Point = self + other (new point, does not modify self)
        """
        return Point(self.x+other.x,self.y+other.y, self.z+other.z)

    def __sub__(self, other):
        """Returns new point = self - other.
        
        Args:
          other: passed point
        
        Returns:
          Point = self - other (new point, does not modify self)
        """
        return Point(self.x-other.x,self.y-other.y, self.z-other.z)

    def __hash__(self):
        """Set the hash of the item to the id number. 
        
        This allows one to make sets of these items.
        """
        return self.id        

    def __mul__(self, factor):
        """Returns new point = self*factor.
        
        Args:
          factor (double): factor to multiply the Point by
        
        Returns:
          Point = self*factor (new point, does not modify self)

        """
        return Point(self.x*factor, self.y*factor, self.z*factor)

    def length(self):
        """Returns the length of this point or vector."""
        res = (self.x**2 + self.y**2)**0.5
        return res

    def make_unit(self):
        """Modifies self vector to make it a unit vector."""
        lenval  = self.length()
        self.x = self.x*1/lenval
        self.y = self.y*1/lenval

    def ang_rad(self):
        """Returns angle in radians, assume point is vector from 0,0."""
        vert = self.x
        horiz = self.y
        radians = atan2(vert, horiz)
        return radians

    def ang_deg(self):
        """Returns angle in degrees, assume point is vector from 0,0."""
        radians = self.ang_rad()
        deg = radians * 180.0 / pi
        return deg

    def rot_ccw_deg(self, ang):
        """Rotates the current vector by ccw degrees about 0,0.
        
        Args:
          ang: angle or rotation in degrees, counter-clockwise (ccw)        
        """
        ang = radians(ang)
        ax = self.y*cos(ang) - self.x*sin(ang)
        rad = self.y*sin(ang) + self.x*cos(ang)
        ax = round(ax, 5)
        rad = round(rad, 5)
        (self.x, self.y) = (rad, ax)

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
        val = 'Point, id %i, (x,y)=(%f,%f)' % (self.id, self.x,self.y)
        return val


class Line(base_classes.Idobj):
    """Stores a line from points p1 to p2.
    
    Arguments:
      p1 (Point): first point
      p2 (Point): second point
      sign (int): 1 or -1, defines direction of line
        +1: p1 -> p2
        -1: p2 -> p1
    
    Attributes:
      _points_ (list of Point) : list of the line's points [p1,p2]
      sign (int): 1 or -1, defines direction of line
        +1: p1 -> p2
        -1: p2 -> p1
      midpt (Point): a mid-point between p1 and p2
      nodes (list of Node): a list of meshed nodes on the line
      faces (list of Face): a list of meshed faces on the line
    """
    
    def __init__(self, p1, p2, sign=1):
        self._points_ = [p1, p2]
        self.sign = sign
        self.midpt = self.mid()
        self.nodes = []
        self.faces = []
        base_classes.Idobj.__init__(self)
        
    def get_name(self):
        """Returns the line name."""
        return 'L'+str(self.id)

    def plot(self, ax):
        """Draws and labels the line onto the passed Matplotlib axis.
        
        Args:
          ax: Matplotlib axis
        """
        lax = [pt.y for pt in self._points_]
        lrad = [pt.x for pt in self._points_]
        ax.plot(lax, lrad)
        pax = self.midpt.y
        prad = self.midpt.x
        ax.annotate(self.get_name(), (pax,prad))

    def length(self):
        """Returns the length of the line."""
        return (self.pt(0)-self.pt(1)).length()

    def set_ediv(self, ediv):
        """Sets the number of element divisions on the line when meshing.
        
        Args:
          ediv (int): number of required elements on this line
        """
        self.ediv = ediv

    def reverse(self):
        """Returns a new line flipped in the opposite direction."""
        a = Line(self._points_[0], self._points_[1], self.sign*-1)
        a.id = self.id
        return a

    def pt(self, ind):
        """Returns the start or end point of the line.
        
        Args:
          ind (int): index of the point we want, 0=start, 1=end
        
        Returns:
          Point: requested Point
        """
        if self.sign == -1:
            return self._points_[not ind]
        else:
            return self._points_[ind]        

    def set_pt(self, ind, pt):
        """Update the line's point at ind to the new pt.
        
        Args:
          ind (int): index of the point we're updating, 0=start, 1=end
          pt (Point): new point we are assigning to Line
        """
        # used to update specific points in the line
        if self.sign == -1:
            self._points_[not ind] = pt
        else:
            self._points_[ind] = pt        
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
        if p1 in other._points_ or p2 in other._points_:
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
        lvect = self.pt(1) - self.pt(0)
        # rotate it ccw by 90 degreeself, and make it a unit
        lvect.rot_ccw_deg(90)
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
          [a,b,c]: list of the above line terms
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
        return [a,b,c]

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
        """Returns an intersecction point on this line of an arc centered at pt
        
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
        [a,b,c] = self.get_abc()
        remainder = a*pt.x + b*pt.y + c
        if remainder < _acc:
            # point is on the line equation, but is it between end pts
            lvect = self.pt(1) - self.pt(0)
            # point0 + lvect*term = pt
            nondim = pt - self.pt(0)
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
        if isinstance(other, Line):
            slope1 = self.get_abc()
            slope2 = other.get_abc()
            if slope1[:2] == slope2[:2]:
                # lines are paralel and will not have an intersection point
                return None
            else:
                # parametric intersection calculation
                # p1 -> p2 = ta parametric, p2 -> p3 tb parametric
                (p1, p2, p3, p4) = (self.pt(0), self.pt(1), other.pt(0), other.pt(1))
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

    def __str__(self):
        """Returns string listing object type, id number, and points"""
        p0 = self.pt(0)
        p1 = self.pt(1)
        val = 'Line, id %i p0: %s p1: %s' % (self.id, p0, p1)
        return val


class Arc(base_classes.Idobj):
    """Makes an arc from points p1 to p2 about center actr.
    
    Arcs should be 90 degrees maximum.
    
    Arguments:
      p1 (Point): first point
      p2 (Point): second point
      actr (Point): arc center
      sign (int): 1 or -1, defines direction of line
        +1: p1 -> p2
        -1: p2 -> p1
    
    Attributes:
      _points_ (list of Point) : list of the line's points [p1,p2]
      sign (int): 1 or -1, defines direction of line
        +1: p1 -> p2
        -1: p2 -> p1
      actr (Point): arc center
      radius (float): radius of the arc
      concavity (str):
        'concave': concave arc
        'convex': convex arc
      midpt (Point): a mid-point between p1 and p2
      nodes (list of Node): a list of meshed nodes on the line
      faces (list of Face): a list of meshed faces on the line
    """

    def __init__(self, p1, p2, actr, sign=1):
        self._points_ = [p1, p2]
        self.sign = sign
        self.actr = actr
        self.radius = (p1-actr).length()
        self.concavity = self.get_concavity()
        self.midpt = self.mid()
        self.nodes = []
        self.faces = []
        base_classes.Idobj.__init__(self)        

    def get_name(self):
        """Returns the arc name."""
        return 'L'+str(self.id)
    
    def plot(self, ax):
        """Draws and labels the arc onto the passed Matplotlib axis.
        
        Args:
          ax: Matplotlib axis
        """
        (pax, prad) = (self.midpt.y, self.midpt.x)
        ax.annotate(self.get_name(), (pax,prad))
        
        # need center, radiuself, ang1, and ang2
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

        # matplotlib assumes ccw arc drawing, calculix assumes cw drawing, reverse order
        a = AArc(xy=[ctr.y, ctr.x], width=2*rad, height=2*rad, angle=0, theta1=ang1, theta2=ang2)
        ax.add_artist(a)

    def length(self):
        """Returns the length of the arc."""
        a = self.pt(0)-self.actr
        b = self.pt(1)-self.actr
        ang = a.ang_bet_rad(b)
        res = self.radius*ang
        return res

    def set_ediv(self, ediv):
        """Sets the number of element divisions on the arc when meshing.
        
        Args:
          ediv (int): number of required elements on this arc
        """
        self.ediv = ediv

    def reverse(self):
        """Returns a new arc flipped in the opposite direction."""
        a = Arc(self._points_[0], self._points_[1], self.sign*-1)
        a.id = self.id
        return a

    def pt(self, ind):
        """Returns the start or end point of the arc.
        
        Args:
          ind (int): index of the point we want, 0=start, 1=end
        
        Returns:
          Point: requested Point
        """
        if self.sign == -1:
            return self._points_[not ind]
        else:
            return self._points_[ind]   
        
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
        if p1 in other._points_ or p2 in other._points_:
            return True
        else:
            return False

    def sgn(self, val):
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
        remainder = (pt.x - self.actr.x)**2 + (pt.y - self.actr.y)**2 - self.radius**2
        if remainder < _acc:
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
                if len(res)==1:
                    return res[0]
                else:
                    return res
        elif isinstance(other, Arc):
            # arc-arc intersection
            # I have not yet written the code here
            pass

    def __str__(self):
        """Returns string listing object type, id number, and points"""
        p0 = self.pt(0)
        p1 = self.pt(1)
        actr = self.actr
        val = 'Arc, id %i p0: %s p1: %s actr: %s' % (self.id, p0, p1, actr)
        return val


class Area(base_classes.Idobj):
    """ Makes an area.
    
    Area is closed in a clockwise direction.
    
    Args:
      p (Part): parent part
      line_list (list): list of Lines and Arcs that close the area
    
    Attributes:
      p (Part): parent part
      closed (boolean): True if closed, False if not
      lines (list): list of lines or arcs that define the area
      matl (Matl): the material fo the area
      etype (str): element type of area. Options are:
        'plstress': plane stress
        'plstrain': plane strain
        'axisym': axisymmetric
      area (double): in-plane part area (ballpark number, arcs not calc right)
      centroid (Point or None): the area centroid point. None before closed.
      nodes (list): child mesh nodes that are in the area
      elements (list): child elements that are in the area
    """
    
    def __init__(self, p, line_list=[]):
        # initialtes the area
        base_classes.Idobj.__init__(self)        
        self.p = p # parent part reference
        self.closed = False
        self.lines = []
        self.matl = None
        self.etype = None
        self.area = 0.0
        self.centroid = None
        self.nodes = []
        self.elements = []
        self.update(line_list)        

    def set_child_ccxtypes(self):
        """Assigns the child element ccx types."""
        if self.etype != None and hasattr(self, 'elements'):
            # if etype is set and elements have been assigned
            for e in self.elements:
                e.set_ccxtype(self.etype)
        
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

    def get_name(self):
        """Returns the area name."""
        return 'A'+str(self.id)
    
    def close(self):
        """Closes the area."""
        # sets closed=True
        # calculates the area center
        self.closed = True
        self.calc_center()        

    def update(self, line_list):
        """Updates the area's line list to the passed list.
        
        The area's line list is replaced with the passed line list.
        If the area is cosed, self.close() will be called.
        
        Args:
          line_list (list): list of Lines or Arcs
        """
        # adds the line list and checks for closure
        self.lines = line_list
        if len(line_list) > 2:
            if self.lines[0].pt(0) == self.lines[-1].pt(1):
                self.close()
                
    def line_from_startpt(self, pt):
        """Returns the line in the area that start at the passed pt.
        
        Args:
          pt (Point): the start point of the line one wants
        
        Returns:
          None or line (Line or Arc): Returns None if no line is found,
            otherwise the correct line or arc which starts with pt will be
            returned.
        """
        # returns the line that starts on the given point
        for line in self.lines:
            if line.pt(0) == pt:
                return line
        return None

    def get_points(self):
        """Returns all the child points in the area.
        
        These are the points that define all lines and arcs in the area.
        Does not return arc center points, only line start and end points.
        Only unique points are returned (no duplicate points).
        Points are generated by iterating through the line_list.
        """
        pts = []
        for line in self.lines:
            p1 = line.pt(0)
            p2 = line.pt(1)
            if p1 not in pts:
                pts.append(p1)
            if p2 not in pts:
                pts.append(p2)
        return pts

    def calc_center(self):
        """Calculates the area centroid, the center of the area.
        
        Attributes area and centroid are set by this method:
          area (double): in-plane part area (ballpark number, arcs not calc right)
          centroid (Point): the area centroid point.
        """
        # get list of points, append first point as last
        pts = self.get_points()
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
        self.area = area
        self.centroid = Point(cx, cy)

    def get_maxlength(self):
        """Returns max distance between points in an area."""
        pts = self.get_points()
        maxlen = 0.0
        # loop through points checking dist to next point
        for ind, p1 in enumerate(pts[:-1]):
            for p2 in pts[ind:]:
                vect = p1 - p2
                dist = vect.length()
                if dist > maxlen:
                    maxlen = dist
        return maxlen

    def add(self, line):
        """Adds line to the area definition.
        
        Line is added to the end of the lines list.
        Area is closed if needed.
        """
        # this adds a line to the area
        self.lines.append(line)
        # check to see if the area is closed
        if line.pt(1) == self.lines[0].pt(0):
            self.close()

    def pt_inside(self, pt):
        """Returns bool telling if pt is inside the area.
        
        Args:
          pt (Point): point we are checking
        
        Returns:
          bool: True if in this area False if not
        """
        # check bounding box first
        # we assume that the point is not on the boundary lines
        pts = self.get_points()
        xs = [p.x for p in pts]
        ys = [p.y for p in pts]
        if (min(xs) <= pt.x <= max(xs)) and (min(ys) <= pt.y <= max(ys)):
            # point is in area bounding box
            vlen = self.get_maxlength()
            pend = Point(pt.x + 0.0001*vlen, pt.y+vlen)
            cutline = Line(pt, pend)
            overlaps = 0
            for line in self.lines:
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
        
    def line_insert(self, lgiven, lnew):
        """Insets line lnew after the given line.
        
        Args:
          lgiven (Line or Arc): the line before the place we will insert a line
          lnew (Line or Arc): the new line we will insert after lgiven
        """
        # this inserts a line after the given line
        if lgiven in self.lines:
            ind = self.lines.index(lgiven)+1
            self.lines.insert(ind,lnew)
            return True
        else:
            return False
    
    def __str__(self):
        """Returns a string listing the area's name."""
        return self.get_name()
        
