from . import base_classes
from . import geometry #point, line, area

#accuracy for small numbers in math below
_acc = .00001

class Part(base_classes.Idobj):
    """This makes a part.

    Args:
      parent: parent FeaModel
    
    Attributes:
      p (FeaModel): parent FeaModel
      cursor (Point): location of a cursor drawing the part
      matl (Matl): the part material
      thickness (float): the part thickness if applicable
        this is only applicable if th epart is plane strain or plane stress
      areas (list): list of Area that make up the part
      left (list): list the parts leftmost lines, they must be vertical
      right (list): list the parts rightmost lines, they must be vertical
      top (list): list the parts top lines, they must be horizontal
      bottom (list): list the parts bottom lines, they must be horizontal
      nodes (list): list of part's nodes
      elements (list): list of part's elements
    """

    def __init__(self, parent):
        self.p = parent
        self.cursor = geometry.Point(0,0)
        self.matl = None
        thickness = None
        self.areas = [] # top area is the buffer
        # make the buffer
        area = self.p.areas.append( geometry.Area(self, []) )
        self.areas.append(area)
        base_classes.Idobj.__init__(self)             
        self.left = None
        self.right = None
        self.top = None
        self.bottom = None
        self.nodes = []
        self.elements= []
        
    def get_item(self, item):
        """"Returns the part's item(s) requested by the passed string.
        
        Args:
          item (str): string requesting item(s)
            Valid examples: 'P0', 'L0', 'left', 'A0'
        
        Returns:
          item(s) or None:
            If items are found they are returned
              If there is only one item it is returned
              If there are multiple items, they are returned as a list
            If no items are found None is returned
        """
        if item in ['left','right','top','bottom']:
            items = getattr(self, item)
            return items
        elif item[0] == 'P':
            # get point
            items = self.get_points()
            num = int(item[1:])
            res = [a for a in items if a.id == num]
            return res[0]
        elif item[0] == 'L':
            # get line
            items = self.get_lines()
            num = int(item[1:])
            items = [a for a in items if a.id == num]
            return items[0]
        elif item[0] == 'A':
            # get area
            items = self.areas
            num = int(item[1:])
            items = [a for a in items if a.id == num]
            return items[0]
        else:
            print('Unknown item! Please pass the name of a point, line or area!')
            return None

    def get_name(self):
        """Returns the part name based on id number."""
        return 'PART'+str(self.id)

    def set_side(self, loc, ind, axis):
        """Sets the part.side to a list of side lines. side=Left,Right,Top,Bottom.
        
        Args:
          loc (string): 'left', 'right', 'top','bottom'
          ind (int): 0 for low, 1 for high
          axis (string): axies to use, 'x' or 'y'
        """
        # loc = 'left', ind = 0, axis = 'y'
        points = self.get_points()
        # sort the points low to high
        points = sorted(points, key=lambda pt: getattr(pt, axis))
        # store the value
        val = getattr(points[ind], axis)
        res = []
        lines = self.get_lines()
        for line in lines:
            if isinstance(line, geometry.Line):
                if getattr(line.pt(0),axis) == val and getattr(line.pt(1),axis) == val:
                    # line is on the left side
                    res.append(line)
        setattr(self, loc, res)
            
    def goto(self, x, y):
        """Moves the part cursor to a location.
        
        If that location has a point at it, use it.
        If not, make a new point at that location.
        
        Args:
          x (float): x-coordinate of the point to go to
          y (float): y-coordinate of the point to go to
        
        Returns:
          self.cursor (Point): returns the updated cursor point
        """
        [pnew, already_exists] = self.make_get_pt(x, y)

        if already_exists:
            if self.areas[-1].closed == True:
                # make a new area if the old area is already closed and we're
                # goint to an existing point
                a = self.p.areas.append( geometry.Area(self, []) )
                self.areas.append( a )
        
        # return cursor
        self.cursor = pnew
        return self.cursor

    def make_get_pt(self, x, y):
        """Gets a point if it exists, makes it if it doesn't. Returns the point.
        
        Use this when you need a point made in the part, and you want to
        use an extant point if one is available.
        
        Args:
          x (float): point x-coordinat
          y (float): point y-coordinate
        
        Returns:
          list:
            list[0]: Point
            list[1]: boolean, True = the point already existed
        """
        pnew = geometry.Point(x,y)
        pts = self.get_points()
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
            pnew = self.p.points.append( pnew )
        return [pnew, pexists]
    
    def make_get_line(self, lnew):
        """Gets a line if it exists, makes it if it doesn't. Returns the line.
        
        Use this when you need a point line in the part, and you want to
        use an extant line if one is available.
        
        Args:
          lnew (Line or Arc): line to make
        
        Returns:
          list:
            list[0]: Line or Arc
            list[1]: boolean, True = the line already existed
        """        
        lnew_pset = set([p.id for p in lnew._points_])
        lines = self.get_lines()
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
            lnew = self.p.lines.append( lnew )        
        return [lnew, lexists]            

    def draw_arc(self, x2, y2, xc, yc):
        """Makes an arc and adds it to the part.

        Current point is the first arc point.
        (x2,y2) is the end point
        (xc,yc) is the arc center
        
        Args:
          x2 (float): arc end point x-coordinate
          y2 (float): arc end point y-coordinate
          xc (float): arc center point x-coordinate
          yc (float): arc center point y-coordinate
        
        Returns:
          list: [arc, arc_start_point, arc_end_point]
        """
        pold = self.cursor
        # make arc center point
        ctr = geometry.Point(xc,yc)
        if ctr in self.p.points:
            ctr = next(x for x in self.p.points if x == ctr)
        else:
            self.p.points.append( ctr )
        # make arc end point
        self.cursor = self.p.points.append( geometry.Point(x2,y2) )
        # make arc
        arc = self.p.lines.append( geometry.Arc(pold, self.cursor, ctr) )
        self.areas[-1].add(arc)
        return [arc, pold, self.cursor]

    def draw_line_delta(self, dx, dy):
        """Draws a line a relative distance, and adds it to the part.

        Args:
          dx (float): x-axis delta distance to draw the line
          dy (float): y-axis delta distance to draw the line
        
        Returns:
          list: [line, point_start, point_end]
        """
        x = self.cursor.x + dx
        y = self.cursor.y + dy
        return self.draw_line_to(x, y)
    
    def draw_line_rad(self, dx_rad):
        """Draws a line a relative radial distance, and adds it to the part.

        Args:
          dx_rad (float): x-axis delta distance to draw the line
        
        Returns:
          list: [line, point_start, point_end]
        """
        return self.draw_line_delta(dx_rad, 0.0)
    
    def draw_line_ax(self, dy_ax):
        """Draws a line a relative axial distance, and adds it to the part.

        Args:
          dy_ax (float): y-axis delta distance to draw the line
        
        Returns:
          list: [line, point_start, point_end]
        """
        return self.draw_line_delta(0.0, dy_ax)    

    def draw_line_to(self, x, y):
        """Draws a line to the given location, and adds it to the part.

        Args:
          x (float): x-axis coordinate of the end point
          y (float): y-axis coordinate of the end point
        
        Returns:
          list: [line, point_start, point_end]
        """
        pold = self.cursor
        [self.cursor, p_already_exists] = self.make_get_pt(x,y)
        [line, l_already_exists] = self.make_get_line( 
            geometry.Line(pold, self.cursor) )
        
        # add the line to the area
        self.areas[-1].add(line)
        
        # check for closure of the area
        if self.areas[-1].closed:
            self.set_side('left',0,'y')
            self.set_side('right',-1,'y')
            self.set_side('top',-1,'x')
            self.set_side('bottom',0,'x')
            self.p.select(self)

        return [line, pold, self.cursor]

    def get_lines(self):
        """Returns a list of all lines in the part."""
        lines = []
        for area in self.areas:
            for line in area.lines:
                if line not in lines:
                    lines.append(line)
        return lines
    
    def get_points(self):
        """Returns a list of all points in the part."""
        points = []
        for area in self.areas:
            for line in area.lines:
                for pt in line._points_:
                    if pt not in points:
                        points.append(pt)
        return points
    
    def get_maxlength(self):
        """Returns the max distance between points in the part."""
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

    def area_from_pt(self, pt):
        """Returns the area that the point is inside.
        
        Args:
          pt (Point): the point we are asking about
        
        Returns:
          Area or None:
            Area is the found area
            None is returned if the point is not in one of this part's areas
        """
        for area in self.areas:
            if area.pt_inside(pt):
                return area
        return None

    def fillet_lines(self, l1, l2, radius):
        """Fillets the given lines in the part.
        
        This inserts an arc in the part tangent to the two given lines.
        
        Args:
          l1 (Line): line that the arc starts on, arc is tangent to the line
          l2 (Line): line that the ar ends on, arc is tangent to the line
          radius (float): arc radius size        
        """
        # this function fillets lines in a part
        # check if the lines are touching
        tmp = self.cursor
        
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
            
            # del old pt, store new points for the arc
            self.p.points.remove(rempt)
            [p1_new,b] = self.make_get_pt( p1_new.x, p1_new.y )
            [ctrpt,b] = self.make_get_pt( ctrpt.x, ctrpt.y )
            [p2_new,b] = self.make_get_pt( p2_new.x, p2_new.y )
            
            # make the new arc
            arc = self.p.lines.append( geometry.Arc(p1_new, p2_new, ctrpt) )
                        
            # edit the adjacent lines to replace the removed pt
            l1.set_pt(1,arc.pt(0))
            l2.set_pt(0,arc.pt(1))
            
            # put the arc in the right location in the area
            for area in self.areas:
                if l1 in area.lines:
                    area.line_insert(l1, arc)
                    print ('Arc inserted into area %i' % (area.id,))
            self.p.sel_children()
                        
        else:
            print ('Cannot fillet! Lines must touch!')
        # reset the cursor to where it should be
        self.cursor = tmp
        
    def point_parents(self, pt):
        """Returns lines that contain pt as their start or end point.
        
        Args:
          pt (Point): the point we are asking about
        
        Returns:
          res (list): list of lines that contain pt
        """
        res = []
        lines = self.get_lines()
        for line in lines:
            if pt in line._points_:
                res.append(line)
        return res
    
    def point_line(self, pt):
        """Returns the line that point pt is on. pt is not a start/end point.
        
        Args:
          pt (Point): the point we are asking about
        
        Returns:
          None or Line or Arc:
            None: if the point is not on one of the lines
            Line or Arc: the line or arc is returned if the pt is on it
        """
        res = None
        lines = self.get_lines()
        for line in lines:
            if line.coincident(pt) == True:
                return line
        return res
    
    def line_parents(self, line):
        """Returns the area parents of the passed line.
        
        Args:
          line (Line or Arc): the line you are asking about
        
        Returns:
          res (list of Area): the areas that contain the line are returned in a
            list
        """
        res = []
        for area in self.areas:
            if line in area.lines:
                res.append(area)
        return res
    
    def cut_line(self, pt, line):
        """Cuts the passed line at the passed point pt.
        
        The passed line is cut into two lines. All areas that included the
        original line are updated.
        
        Args:
          line (Line): the line to cut
          pt (Point): the location on the line that we will cut it
        
        Returns:
          list: [pnew, lnew]
            pnew: the new point we created to cut the original line
            lnew: the new line we created, the end half of the orignal line
        """
        p = geometry.Point(pt.x,pt.y)
        pnew = None
        if p in self.p.points:
            pnew = [pt for pt in self.p.points if pt == p][0]
        else:
            pnew = self.p.points.append( p )
        # below is code to split line, requires defined point and 
        pend = line.pt(1)
        line.set_pt(1, pnew)
        lnew = self.p.lines.append( geometry.Line(pnew, pend) ) 
        areas = self.line_parents(line)
        # insert the new line into existing areas
        for area in areas:
            area.line_insert(line, lnew)
        return [pnew, lnew]

    def cut(self, pt, cvect):
        """Cuts the part at point pt using cvect as a cutting vector.

        This function doesn't return anything.
        Instead, it updates all part areas, lines etc after cutting.
        So if you cut one area into two pieces, after calling this function
        the part will have two areas.
        
        Args:
          pt (Point): the locatin we are cutting from
          cvect (Point): the vector direction of the cut from pt
        """
        # this cuts the part, through one or more areas
        cvect.make_unit()
        
        # is the point an end point of a line?
        lines_excluded = self.point_parents(pt)
        
        # if the pt is on an existing line, add that line to the excl list
        if len(lines_excluded) == 0:
            L1 = self.point_line(pt)
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
            [pnew, lnew] = self.cut_line(pt, line_origin)
            # set the first point of our cut line to be this actual point
            linepts[0] = pnew
        
        # make the tool cut line
        vsize = self.get_maxlength()
        endpt = pt + cvect*vsize
        cutline = geometry.Line(pt,endpt)
        
        # find the intersections
        overlaps = []
        lines = self.get_lines()
        for l in lines:
            if l not in lines_excluded:
                if isinstance(l, geometry.Line):
                    newpt = l.intersects(cutline)
                    if isinstance(newpt, geometry.Point):
                        areas = len(self.line_parents(l))
                        
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
                [pnew, lnew] = self.cut_line(cut['pt'], cut['line'])
            # adds the intersection point to the list of cut line points
            linepts.append(pnew)
            p1 = linepts[-2]
            p2 = linepts[-1]
            pav = p1+p2
            pav = pav*0.5
            
            # get area to cut
            area = self.area_from_pt(pav)
            # make the cut line
            cut_line = self.p.lines.append( geometry.Line(p1, p2) )
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
            anew = self.p.areas.append( geometry.Area(self, alist_other) )
            self.areas.append(anew)

            if cut['areas'] == 1:
                # break out of for loop, we hit a part edge
                break
        
    def chunk_area(self, area):
        """Cuts the passed area into regular smaller areas.
        
        The cgx mesher only accepts areas which are 3-5 sides
        so one may need to call this before using that mesher.
        Cuts are made perpendicular to tangent points or at
        internal corners.
        At internal corners two perpendicular cuts are made.
        
        Args:
          area (Area): the area to cut into smaller areas
        """
        # store the cuts first, then cut after
        cuts = [] # each item is a dict with a pt and vect in it
        for ind, line in enumerate(area.lines):
            L1 = area.lines[ind-1]
            L2 = line
            pt = L1.pt(1)
            v1 = L1.get_perp_vec(pt)
            v2 = L2.get_perp_vec(pt)
            # flip these vectors later to make them cut the area(s)
            ang = v1.ang_bet_deg(v2)
            d = {}
            if ang == 0.0:
                # tangent
                d = {'pt':pt, 'vect':v1*-1}
                cuts.append(d)
            elif ang > 0:
                # internal corner
                d = {'pt':pt, 'vect':v1*-1}
                cuts.append(d)
                d = {'pt':pt, 'vect':v2*-1}
                cuts.append(d)
            elif ang < 0:
                # external corner
                # do not split these
                pass
        
        # do the cuts
        for cut in cuts:
            print('--------------------')
            print('Cut pt:',cut['pt'])
            print('Cut vect: ',cut['vect'])
            self.cut(cut['pt'], cut['vect'])
            
    def chunk(self):
        """Chunks all areas in the part."""
        for area in self.areas:
            if area.closed:
                if len(area.lines) > 5:
                    # we need to chunk this area, it has too many lines to mesh
                    self.chunk_area(area)
                else:
                    print('Area %i was not chunked because it had < 6 sides.'%(area.id))
        # store the left, right, top, and bottom lines        
        self.set_side('left',0,'y')
        self.set_side('right',-1,'y')
        self.set_side('top',-1,'x')
        self.set_side('bottom',0,'x')
        self.p.select(self)
        
    def set_matl(self, mat):
        """Sets the part material.
        
        Args:
          matl (Matl): the material assigned to the part
        """
        self.matl = mat
            
    def __str__(self):
        """Returns string listing object type, id number and name."""
        val = 'Part, id=%i name=%s' % (self.id, self.get_name())
        return val
