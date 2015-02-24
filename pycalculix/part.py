"""This module stores the Part class. It is used to make 2D parts.
"""

from . import base_classes
from . import geometry #point, line, area

#accuracy for small numbers in math below
ACC = .00001

class Part(base_classes.Idobj):
    """This makes a part.

    Args:
        parent: parent FeaModel

    Attributes:
        __fea (FeaModel): parent FeaModel
        points (list): list or part points, excludes arc centers
        allpoints (list): list or part points, includes arc centers
        lines (list): list of all Line and Arc that make the part
        signlines (list): list of all SignLine and SignArc that make the part
        __cursor (Point): location of a cursor drawing the part
        __holemode (bool): if True, lines will be added to holes, otherwise,
            they'll be added to areas
        areas (list): list of Area that make up the part
        left (list): list the parts leftmost lines, they must be vertical
        right (list): list the parts rightmost lines, they must be vertical
        top (list): list the parts top lines, they must be horizontal
        bottom (list): list the parts bottom lines, they must be horizontal
        nodes (list): list of part's nodes
        elements (list): list of part's elements
    """

    def __init__(self, feamodel):
        self.__fea = feamodel
        self.__cursor = geometry.Point(0, 0)
        self.areas = [] # top area is the buffer
        # make the buffer
        area = self.__fea.areas.append(geometry.Area(self, []))
        self.areas.append(area)
        base_classes.Idobj.__init__(self)
        self.left = None
        self.right = None
        self.top = None
        self.bottom = None
        self.nodes = []
        self.elements = []
        self.__holemode = False
        self.__fea.parts.append(self)

    def __hash__(self):
        """Returns the item's id as its hash."""
        return self.id

    @property
    def lines(self):
        """Returns list of part lines."""
        lines = set()
        for area in self.areas:
            lines.update(area.lines)
        return list(lines)

    @property
    def signlines(self):
        """Returns list of part signline and signarc."""
        lines = set()
        for area in self.areas:
            lines.update(area.signlines)
        return list(lines)

    @property
    def points(self):
        """Returns list of part points, excludes arc centers."""
        points = set()
        lines = self.lines
        for line in lines:
            points.update(line.points)
        return list(points)

    @property
    def allpoints(self):
        """Returns list of part points, includes arc centers."""
        points = set()
        lines = self.lines
        for line in lines:
            points.update(line.allpoints)
        return list(points)

    def get_item(self, item):
        """"Returns the part's item(s) requested by the passed string.

        Args:
            item (str): string requesting item(s)

                * Valid examples: 'P0', 'L0', 'left', 'A0'

        Returns:
            item(s) or None: If items are found they are returned

                * If there is only one item it is returned
                * If there are multiple items, they are returned as a list
                * If no items are found None is returned
        """
        if item in ['left', 'right', 'top', 'bottom']:
            items = getattr(self, item)
            return items
        elif item[0] == 'P':
            # get point
            items = self.points
            num = int(item[1:])
            res = [a for a in items if a.id == num]
            return res[0]
        elif item[0] == 'L':
            # get line
            items = self.signlines
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

    def __set_side(self, side):
        """Sets the part.side to a list of lines on that side of the part.

        Used to set the part.left, part.right, part.top, part.bottom sides.

        Args:
            side (string): 'left', 'right', 'top','bottom'
        """
        # set index and axis, ind=0 is low side, ind=-1 is high side
        inds = {'left':0, 'right':-1, 'top':-1, 'bottom':0}
        axes = {'left':'y', 'right':'y', 'top':'x', 'bottom':'x'}
        ind = inds[side]
        axis = axes[side]

        # loc = 'left', ind = 0, axis = 'y'
        points = self.points
        # sort the points low to high
        points = sorted(points, key=lambda pt: getattr(pt, axis))
        # store the value
        val = getattr(points[ind], axis)
        res = []
        lines = self.signlines
        for sline in lines:
            if isinstance(sline, geometry.SignLine):
                vals = [getattr(pt, axis) for pt in sline.points]
                if vals == [val, val]:
                    # line is on the left side
                    res.append(sline)
        setattr(self, side, res)

    def goto(self, x, y, holemode=False):
        """Moves the part cursor to a location.

        If that location has a point at it, use it.
        If not, make a new point at that location.

        Args:
            x (float): x-coordinate of the point to go to
            y (float): y-coordinate of the point to go to
            holemode (bool): if True, we start drawing a hole here, otherwise
                we start drawing an area

        Returns:
            self.__cursor (Point): returns the updated cursor point
        """
        [pnew, already_exists] = self.__make_get_pt(x, y)

        if already_exists:
            if self.areas[-1].closed == True:
                # make a new area if the old area is already closed and we're
                # going to an existing point
                area = self.__fea.areas.append(geometry.Area(self, []))
                self.areas.append(area)

        # return cursor
        self.__cursor = pnew
        self.__holemode = holemode
        return self.__cursor

    def __get_point(self, point):
        """Returns point if found, None otherwise."""
        points = self.allpoints
        found_point = None
        for apoint in points:
            dist = point - apoint
            dist = dist.length()
            if dist < ACC:
                # point already exists in part, use it
                found_point = apoint
                break
        return found_point

    def __make_get_pt(self, x, y):
        """Gets a point if it exists, makes it if it doesn't. Returns the point.

        Use this when you need a point made in the part, and you want to
        use an extant point if one is available.

        Args:
            x (float): point x-coordinate
            y (float): point y-coordinate

        Returns:
            list:
                list[0]: Point
                list[1]: boolean, True = the point already existed
        """
        thept = geometry.Point(x, y)
        pfound = self.__get_point(thept)
        pexists = True
        if pfound == None:
            pfound = self.__fea.points.append(thept)
            pexists = False
        return [pfound, pexists]

    def __make_get_sline(self, lnew):
        """Returns a signed line or arc, makes it if it needs to.

        Args:
            lnew (Line or Arc or SignLine or SignArc): Line or Arc to make

        Returns:
            list:
                list[0]: SignLine or SignArc
                list[1]: boolean, True = the line already existed
        """
        lpos = lnew.signed_copy(1)
        lneg = lnew.signed_copy(-1)
        # get part's signed lines
        lines = self.signlines
        lexists = False
        for sline in lines:
            if lpos == sline:
                lexists = True
                signline_new = sline
                break
            elif lneg == sline:
                lexists = True
                signline_new = sline.signed_copy(-1)
                # set the edge property to False, we now have two slines here
                signline_new.edge = False
                signline_new = self.__fea.signlines.append(signline_new)
                break
        else:
            # fired when we haven't broken out of the loop, the line is new
            saved_line = self.__fea.lines.append(lnew)
            saved_line.save_to_points()
            signline_new = saved_line.signed_copy(1)
            signline_new = self.__fea.signlines.append(signline_new)
        signline_new.set_parent(self.areas[-1])
        signline_new.line.add_signline(signline_new)
        return [signline_new, lexists]

    def draw_hole(self, center_x, center_y, radius, num_arcs=4, filled=False):
        """Makes a hole in the part.

        Args:
            center_x (float): x-axis hole center
            center_y (float): y-axis hole center
            radius (float): hole radius
            num_arcs (int): number of arcs to use
            filled (bool): whether to fill the hole

                * True: makes a new area in the part

        Returns:
            hole_lines (list or None): list of hole SignLine or SignArc

                * Returns None if hole was not made.
        """
        center = geometry.Point(center_x, center_y)
        area = self.__area_from_pt(center)
        if area == None:
            print("You can't make a hole here until there's an area here!")
            return None
        else:
            rvect = geometry.Point(0, radius)
            pold = center + rvect
            self.goto(pold.x, pold.y, holemode=True)
            for ind in range(num_arcs):
                rvect.rot_ccw_deg(360/num_arcs)
                pnew = center + rvect
                self.draw_arc(pnew.x, pnew.y, center.x, center.y)
                pold = pnew
            # make new area
            if filled:
                lines = list(area.holes[-1])
                lines.reverse() #reverse the order
                lines = [L.signed_copy(-1) for L in lines] #reverse orientation
                # store them in the feamodel
                lines = [self.__make_get_sline(line)[0] for line in lines]
                anew = self.__fea.areas.append(geometry.Area(self, lines))
                self.areas.append(anew)
            return area.holes[-1]

    def draw_arc_angle(self, degrees_ccw, center_x, center_y):
        """Makes an arc and adds it to the part.

        | Current point is the first arc point.
        | degrees_ccw is the swept angle in degrees, counterclockwise
        | (center_x, center_y) is the arc center

        | Degrees: Traversed angle of arc must be < 180 degrees

        Args:
            degrees_ccw (float): arc swept angle in degrees, counterclockwise
            center_x (float): arc center point x-coordinate
            center_y (float): arc center point y-coordinate

        Returns:
            list: [arc, arc_start_point, arc_end_point]
        """
        center = geometry.Point(center_x, center_y)
        radius_vector = self.__cursor - center
        radius_vector.rot_ccw_deg(degrees_ccw)
        end = center + radius_vector
        return self.draw_arc(end.x, end.y, center_x, center_y)

    def draw_arc(self, end_x, end_y, center_x, center_y):
        """Makes an arc and adds it to the part.

        | Current point is the first arc point.
        | (end_x, end_y) is the end point
        | (center_x, center_y) is the arc center

        | Degrees: Traversed angle of arc must be < 180 degrees
        | Radians: Traversed angle of arc must be < Pi

        Args:
            end_x (float): arc end point x-coordinate
            end_y (float): arc end point y-coordinate
            center_x (float): arc center point x-coordinate
            center_y (float): arc center point y-coordinate

        Returns:
            list: [arc, arc_start_point, arc_end_point]
        """
        pold = self.__cursor
        # make arc center point
        ctr = self.__make_get_pt(center_x, center_y)[0]
        # make arc end point
        self.__cursor = self.__make_get_pt(end_x, end_y)[0]

        # make arc
        arc = self.__make_get_sline(geometry.Arc(pold, self.__cursor, ctr))[0]

        if self.__holemode:
            area = self.__area_from_pt(self.__cursor)
            if area != None:
                closed = area.add_hole_sline(arc)
                if closed:
                    self.__holemode = False
            else:
                print('You must have a closed area here before making a hole!')
        else:
            self.areas[-1].add_sline(arc)
        return [arc, pold, self.__cursor]

    def draw_line_delta(self, delta_x, delta_y):
        """Draws a line a relative distance, and adds it to the part.

        Args:
            delta_x (float): x-axis delta distance to draw the line
            delta_y (float): y-axis delta distance to draw the line

        Returns:
            list: [line, point_start, point_end]
        """
        x = self.__cursor.x + delta_x
        y = self.__cursor.y + delta_y
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
            list: [SignLine, point_start, point_end]
        """
        pold = self.__cursor
        self.__cursor = self.__make_get_pt(x, y)[0]
        sline = self.__make_get_sline(geometry.Line(pold, self.__cursor))[0]

        if self.__holemode:
            area = self.__area_from_pt(self.__cursor)
            if area != None:
                closed = area.add_hole_sline(sline)
                if closed:
                    self.__holemode = False
            else:
                print('You must have a closed area here before making a hole!')
        else:
            # drawing in last area
            self.areas[-1].add_sline(sline)
            # check for closure of the area
            if self.areas[-1].closed:
                self.__set_side('left')
                self.__set_side('right')
                self.__set_side('top')
                self.__set_side('bottom')
        return [sline, pold, self.__cursor]

    def __get_maxlength(self):
        """Returns the max distance between points in the part."""
        points = self.points
        maxlen = 0.0
        # loop through points checking dist to next point
        for ind, point_1 in enumerate(points[:-1]):
            for point_2 in points[ind:]:
                vect = point_1 - point_2
                dist = vect.length()
                if dist > maxlen:
                    maxlen = dist
        return maxlen

    def __area_from_pt(self, point):
        """Returns the area that the point is inside.

        Args:
            point (Point): the point we are asking about

        Returns:
            Area or None:
                Area is the found area
                None is returned if the point is not in one of this part's areas
        """
        for area in self.areas:
            if area.contains_point(point):
                return area
        return None

    def fillet_lines(self, line1, line2, radius):
        """Fillets the given lines in the part.

        This inserts an arc in the part tangent to the two given lines.

        Args:
            line1 (Line): line that the arc starts on, arc is tangent to the line
            line2 (Line): line that the ar ends on, arc is tangent to the line
            radius (float): arc radius size
        
        Returns:
            list: [arc, start_point, end_point]
        """
        # this function fillets lines in a part
        # check if the lines are touching
        tmp = self.__cursor

        if line1.touches(line2):
            # offset the lines, assuming area is being traced clockwise
            # get the intersection point
            magnitude = radius
            l1_off = line1.offset(magnitude)
            l2_off = line2.offset(magnitude)
            ctrpt = l1_off.intersects(l2_off)
            if ctrpt == None:
                # flip the offset direction if lines don't intersect
                magnitude = -radius
                l1_off = line1.offset(magnitude)
                l2_off = line2.offset(magnitude)
                ctrpt = l1_off.intersects(l2_off)

            # now we have an intersecting point
            print('Arc center pt is: ', ctrpt)
            p1_new = line1.arc_tang_intersection(ctrpt, magnitude)
            p2_new = line2.arc_tang_intersection(ctrpt, magnitude)
            rempt = line1.pt(1)

            p1_new = self.__make_get_pt(p1_new.x, p1_new.y)[0]
            ctrpt = self.__make_get_pt(ctrpt.x, ctrpt.y)[0]
            p2_new = self.__make_get_pt(p2_new.x, p2_new.y)[0]

            # make the new arc
            arc = self.__make_get_sline(geometry.Arc(p1_new, p2_new, ctrpt))[0]

            # edit the adjacent lines to replace the removed pt
            line1.set_pt(1, arc.pt(0))
            line2.set_pt(0, arc.pt(1))
            # del old pt, store new points for the arc
            self.__fea.points.remove(rempt)

            # put the arc in the right location in the area
            for area in self.areas:
                if line1 in area.signlines:
                    area.line_insert(line1, arc)
                    print('Arc inserted into area %i' % (area.id))
                    return [arc, arc.pt(0), arc.pt(1)]
        else:
            print('Cannot fillet! Lines must touch!')
        # reset the cursor to where it should be
        self.__cursor = tmp

    def __cut_line(self, point, line):
        """Cuts the passed line at the passed point.

        The passed line is cut into two lines. All areas that included the
        original line are updated.

        Args:
            line (Line): the line to cut
            point (Point): the location on the line that we will cut it

        Returns:
            list: [pnew, lnew]
                pnew: the new point we created to cut the original line
                lnew: the new line we created, the end half of the orignal line
        """
        pnew = self.__make_get_pt(point.x, point.y)[0]
        pend = line.pt(1)
        line.set_pt(1, pnew)
        lnew = self.__make_get_sline(geometry.Line(pnew, pend))[0]
        # insert the new line into existing areas
        areas = line.areas
        for area in areas:
            area_line = [lin for lin in area.signlines if lin.line == line]
            area_line = area_line[0]
            if area_line.sign == 1:
                # cutting line in clockwise area, where line is pos
                area.line_insert(area_line, lnew)
            elif area_line.sign == -1:
                # cutting line in clockwise area, where line is neg
                lnew_rev = self.__make_get_sline(lnew.signed_copy(-1))[0]
                area.line_insert(area_line, lnew_rev, after=False)
        return [pnew, lnew]

    def __cut_area(self, area, start_pt, end_pt):
        """Cuts the part area from start_pt to end_pt."""
        # select the line portions that define the areas
        # list[:low] excludes low index
        # list [high:] includes high index
        #    we want the line which start with the point
        lpre_start = area.line_from_startpt(start_pt)
        lpre_end = area.line_from_startpt(end_pt)
        istart = area.exlines.index(lpre_start)
        iend = area.exlines.index(lpre_end)
        low = min(istart, iend)
        high = max(istart, iend)

        # lists of lines for areas
        beg = area.exlines[:low]
        mid = area.exlines[low:high]
        end = area.exlines[high:]

        # make cut line for [beg + cut + end] area
        start_pt = mid[0].pt(0)
        end_pt = mid[-1].pt(1)
        fwd = geometry.Line(start_pt, end_pt)
        rev = geometry.Line(end_pt, start_pt)

        # update existing area
        cline = self.__make_get_sline(fwd)[0]
        alist_curr = beg + [cline] + end
        area.update(alist_curr)

        # make new area
        cline_rev = self.__make_get_sline(rev)[0]
        alist_other = mid + [cline_rev]
        anew = self.__fea.areas.append(geometry.Area(self, alist_other))
        self.areas.append(anew)

    def __merge_hole(self, area, start_pt, end_pt):
        """Merges the hole at start_pt into the area."""
        lpre_start = area.line_from_startpt(start_pt)
        hole_line = area.line_from_startpt(end_pt)
        ind = area.exlines.index(lpre_start)

        # store sections of the area
        beg = area.exlines[:ind]
        end = area.exlines[ind:]
        mid = []
        for hole in area.holes:
            for sline in hole:
                if sline == hole_line:
                    mid = hole
                    break
            if mid != []:
                break

        fwd = geometry.Line(start_pt, end_pt)
        fwd_sline = self.__make_get_sline(fwd)[0]
        rev_sline = fwd_sline.signed_copy(-1)
        rev_sline = self.__fea.signlines.append(rev_sline)
        rev_sline.line.add_signline(rev_sline)
        alist_curr = beg + [fwd_sline] + mid + [rev_sline] + end
        area.holes.remove(mid)
        area.update(alist_curr)
        print('Hole merged into area %s' % area)

    def __cut_with_line(self, cutline):
        """Cuts the part using the passed line."""
        # find all intersections
        lines = self.lines
        points = set()
        # add origin point if it is in the part
        if cutline.pt(0) in self.points:
            points.add(geometry.Point(cutline.pt(0).x, cutline.pt(0).y))
        for line in lines:
            newpt = line.intersects(cutline)
            if newpt != None:
                points.add(newpt)

        # loop through intersection points, storing distance and lines to cut
        points = list(points)
        for (ind, point) in enumerate(points):
            dist = point - cutline.pt(0)
            dist = dist.length()
            pdict = {'dist': dist}
            realpt = self.__get_point(point)
            if realpt == None:
                # point does not exist, we will need to cut a line
                realpt = point
                for line in lines:
                    if line.coincident(realpt):
                        pdict['line'] = line
                        break
            pdict['point'] = realpt
            points[ind] = pdict

        # sort the points by dist, lowest to highest
        points = sorted(points, key=lambda k: k['dist'])

        # loop through the points cutting areas
        for ind in range(len(points)):
            pdict = points[ind]
            start_pt = points[ind]['point']
            print(start_pt)
            if 'line' in pdict:
                # cut the line and point to the real new point
                pnew = self.__cut_line(start_pt, pdict['line'])[0]
                points[ind]['point'] = pnew
                start_pt = pnew
            end_pt = None
            pavg = None
            area = None
            if ind > 0:
                # find the area we're working on
                end_pt = points[ind-1]['point']
                pavg = start_pt + end_pt
                pavg = pavg*0.5
                area = self.__area_from_pt(pavg)

                start_hole = start_pt in area.holepoints
                end_hole = end_pt in area.holepoints
                if start_hole and end_hole and area == None:
                    # stop cutting if we are trying to cut through a holes
                    break
                elif start_hole and end_hole and area != None:
                    break
                    # stop cutting if we are trying to join holes
                if end_hole == True:
                    self.__merge_hole(area, start_pt, end_pt)
                else:
                    self.__cut_area(area, start_pt, end_pt)

    def __cut_with_vect(self, point, cvect):
        """Cuts the part at the point using cvect as a cutting vector.

        This function doesn't return anything.
        Instead, it updates all part areas, lines etc after cutting.
        So if you cut one area into two pieces, after calling this function
        the part will have two areas.

        Args:
            point (Point): the location we are cutting from
            cvect (Point): the vector direction of the cut from pt
        """
        cvect.make_unit()
        vsize = self.__get_maxlength()
        endpt = point + cvect*vsize
        cutline = geometry.Line(point, endpt)
        self.__cut_with_line(cutline)

    def __chunk_area(self, area, mode):
        """Cuts the passed area into regular smaller areas.

        The cgx mesher only accepts areas which are 3-5 sides
        so one may need to call this before using that mesher.
        Cuts are made perpendicular to tangent points or at
        internal corners.
        At internal corners two perpendicular cuts are made.

        Args:
            area (Area): the area to cut into smaller areas
            mode (str): 'both', 'holes' or 'ext' chunks the area using the
                points form this set. See part.chunk
        """
        # store the cuts first, then cut after
        cuts = [] # each item is a dict with a pt and vect in it
        loops = []
        if mode == 'holes':
            loops = area.holes
        elif mode == 'ext':
            loops = [area.exlines]
        elif mode == 'both':
            loops = area.holes + [area.exlines]
        for loop in loops:
            for ind, line in enumerate(loop):
                line_pre = loop[ind-1]
                line_post = line
                point = line_pre.pt(1)
                perp1 = line_pre.get_perp_vec(point)
                perp2 = line_post.get_perp_vec(point)
                # flip these vectors later to make them cut the area(s)
                ang = perp1.ang_bet_deg(perp2)
                cut = {}
                if ang == 0.0:
                    # tangent
                    cut = {'pt':point, 'vect':perp1*-1}
                    cuts.append(cut)
                elif ang > 0:
                    # internal corner
                    cut = {'pt':point, 'vect':perp1*-1}
                    cuts.append(cut)
                    cut = {'pt':point, 'vect':perp2*-1}
                    cuts.append(cut)
                elif ang < 0:
                    # external corner
                    # do not split these
                    pass

        # do the cuts
        for cut in cuts:
            print('--------------------')
            print('Cut pt:', cut['pt'])
            print('Cut vect: ', cut['vect'])
            self.__cut_with_vect(cut['pt'], cut['vect'])

    def chunk(self, mode='both'):
        """Chunks all areas in the part.

        Args:
            mode (str): area chunking mode

                - 'both': cuts areas using holes and exterior points
                - 'holes': cut areas using holes points only
                - 'ext': cut areas using exterior points only
        """
        for area in self.areas:
            if area.closed:
                min_sides = 5
                has_holes = len(area.holes) > 0
                ext_gr = len(area.exlines) > min_sides
                both_false = (has_holes == False and ext_gr == False)
                if mode == 'holes' and has_holes:
                    self.__chunk_area(area, mode)
                elif mode == 'both' and (has_holes or ext_gr):
                    self.__chunk_area(area, mode)
                elif mode == 'ext' and ext_gr:
                    self.__chunk_area(area, mode)
                else:
                    aname = area.get_name()
                    val = 'Area %s was not chunked because it had' % aname
                    adder = ''
                    if mode == 'both' and both_false:
                        adder = '<= %i lines and no holes.' % min_sides
                    elif has_holes == False and (mode in ['both', 'holes']):
                        adder = 'no holes.'
                    elif ext_gr == False and (mode in ['both', 'ext']):
                        adder = '<= %i lines.' % min_sides
                    print('%s %s' % (val, adder))

        # store the left, right, top, and bottom lines
        self.__set_side('left')
        self.__set_side('right')
        self.__set_side('top')
        self.__set_side('bottom')

    def __str__(self):
        """Returns string listing object type, id number and name."""
        val = 'Part, id=%i name=%s' % (self.id, self.get_name())
        return val
