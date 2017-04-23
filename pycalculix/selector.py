"""This module stores the View and Selector objects

Selector is used in FeaModel.focus to store the model's selected sets.
View is used in FeaModel.view to plot items and set orientation.
"""

from . import partmodule
from . import geometry
from . import mesh
from . import components

from matplotlib.collections import PatchCollection  # element plotting
from matplotlib.patches import Polygon  # needed for plotting elements
from matplotlib.patches import Arc as AArc
from numpy.linalg import det # needed for arc plotting

# element colors, 0-1, 0=black, 1=whate
ECOLOR = '.4'
FCOLOR = '.9'
CMAP = 'jet' # for results
GEOM_CMAP = 'Pastel1' # color map for parts or areas

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


class View(object):
    """Makes a view object which allows plotting of all entities.

    View orientation is also stored.

    Args:
        orientation (str): view orientation
        focus (Selector): the model's selected items

    Attributes:
        orientation (str): view orientation
    """

    def __init__(self, focus):
        self.__focus = focus
        self.set_orientation(focus.orientation)

    @staticmethod
    def __list_x(items):
        """Private function to return list of xs"""
        return [item.x for item in items]

    @staticmethod
    def __list_y(items):
        """Private function to return list of ys"""
        return [item.y for item in items]

    @staticmethod
    def __coords_xy(items):
        """Private function to return list of [[x1,y1],[x2,y2]...]"""
        return [[item.x, item.y] for item in items]

    @staticmethod
    def __coords_yx(items):
        """Private function to return list of [[y1,x1],[y2,x2]...]"""
        return [[item.y, item.x] for item in items]

    @staticmethod
    def __coord_xy(item):
        """Private function to return Point with x horiz y vert"""
        return geometry.Point(item.x,item.y,item.z)

    @staticmethod
    def __coord_yx(item):
        """Private function to return Point with x horiz y vert"""
        return geometry.Point(item.y,item.x,item.z)

    def __make_plot_arrows(self, figure, loads):
        """Plots the passed loads on to the passed figure using set orientation

        Args:
            figure: matplotlib figure to plot on
            loads (list of load): pressure, displacement loads

                All passed loads must be of the same type
        """
        pass

    def plot_nodes(self, axis, label=True):
        """Plots selected nodes on the passed axis
        
        Args:
            axis: matplotlib axis to plot on
            label (bool): True to plot labels, False otherwise
        
        Returns:
            list: [horizontals, verticals] used to set plot bounds later
        """
        items = self.__focus.nodes
        horiz, vert = self.__horiz(items), self.__vert(items)
        axis.scatter(horiz, vert, s=7, color='black')
        if label:
            names = [node.get_name() for node in items]
            nlist = list(zip(names, horiz, vert))
            for (name, hval, vval) in nlist:
                axis.annotate(name, (hval, vval))
        return [horiz, vert]

    def plot_elements(self, axis, label=True):
        """Plots selected elements on the passed axis
        
        Args:
            axis: matplotlib axis to plot on
            label (bool): True to plot labels, False otherwise
        
        Returns:
            list: [horizontals, verticals] used to set plot bounds later
        """
        items = self.__focus.elements
        polys = []
        bound_nodes = set()
        label_nodes = self.__focus.nodes
        for element in items:
            enodes = element.nodes
            bound_nodes.update(enodes)
            corner_nodes = [node for node in enodes if node.order == 1]
            coords = self.__coords(corner_nodes)
            poly = Polygon(coords, closed=True)
            polys.append(poly)
        coll = PatchCollection(polys, facecolors=FCOLOR, edgecolors=ECOLOR)
        axis.add_collection(coll)

        # plot element numbers
        if label:
            for element in items:
                axis.text(element.center.y, element.center.x,
                          element.get_name(), ha='center', va='center')

        # return horiz, vert to set bounds with
        horiz, vert = self.__horiz(bound_nodes), self.__vert(bound_nodes)
        return [horiz, vert]

    def plot_points(self, axis, label=True):
        """Plots selected points on the passed axis
        
        Args:
            axis: matplotlib axis to plot on
            label (bool): True to plot labels, False otherwise
        
        Returns:
            list: [horizontals, verticals] used to set plot bounds later
        """
        items = self.__focus.points
        horiz, vert = self.__horiz(items), self.__vert(items)
        axis.scatter(horiz, vert)
        if label:
            names = [point.get_name() for point in items]
            nlist = list(zip(names, horiz, vert))
            for (name, hval, vval) in nlist:
                axis.annotate(name, (hval, vval))
        return [horiz, vert]

    def plot_lines(self, axis, label=True):
        """Plots selected lines on the passed axis
        
        Args:
            axis: matplotlib axis to plot on
            label (bool): True to plot labels, False otherwise
        
        Returns:
            list: [horizontals, verticals] used to set plot bounds later
        """
        items = self.__focus.lines
        points = set()
        for item in items:
            if isinstance(item, geometry.SignLine):
                points.update(item.points)
                horiz = self.__horiz(item.points)
                vert = self.__vert(item.points)
                axis.plot(horiz, vert)
            elif isinstance(item, geometry.SignArc):
                points.update(item.points)
                ctr = self.__coord(item.actr)
                rad = item.radius
                pt1 = self.__coord(item.pt(0))
                pt2 = self.__coord(item.pt(1))
                pstart, pend = [pt1, pt2]
                det_val = det([[pt1.x, pt1.y], [pt2.x, pt2.y]])
                # pos is ccw, neg is cw
                if det_val < 0:
                    pstart, pend = [pt2, pt1]
                ang1, ang2 = [pstart - ctr, pend - ctr]
                ang1, ang2 = [ang1.ang_deg(), ang2.ang_deg()]            
                # matplotlib assumes ccw arc drawing, calculix can be cw or ccw
                artist = AArc(xy=[ctr.x, ctr.y], width=2*rad, height=2*rad,
                              angle=0, theta1=ang1, theta2=ang2)
                axis.add_artist(artist)

        # label the lines
        if label:
            for item in items:
                coord = self.__coord(item.midpt)
                lvect = item.get_perp_vec(coord)
                lvect = self.__coord(lvect)
                ang = lvect.ang_deg()
                h_align, v_align = get_text_hv(ang)
                axis.text(coord.x, coord.y, item.get_name(),
                          ha=h_align, va=v_align)
        #return points
        horiz, vert = self.__horiz(points), self.__vert(points)
        return [horiz, vert]

    def plot_areas(self, axis, label=True):
        pass

    def plot_parts(self, axis, label=True):
        pass

    @staticmethod
    def set_bounds(plt, axials, radials):
        """Sets the axial and radial bounds of the shown plot."""
        vert = max(radials) - min(radials)
        horiz = max(axials) - min(axials)
        vadder = (vert)/5
        hadder = (horiz)/5
        (vmax, vmin) = (max(radials)+vadder, min(radials)-vadder)
        (hmax, hmin) = (max(axials)+hadder, min(axials)-hadder)
        plt.xlim(hmin, hmax)
        plt.ylim(vmin, vmax)

    def set_orientation(self, orientation):
        """Sets the view orientation to 'xy' or 'yx' order is horiz, vert

        Args:
            orientation (str): the view orientation

                'xy': default, x is horizontal, y is vertical
                'yx': y is horizontal, x is vertical

        """
        self.orientation = orientation
        print('Orientation set to '+orientation)
        if self.orientation == 'xy':
            self.__horiz = self.__list_x
            self.__vert = self.__list_y
            self.__coords = self.__coords_xy
            self.__coord = self.__coord_xy
        elif self.orientation == 'xy':
            self.__horiz = self.__list_y
            self.__vert = self.__list_x
            self.__coords = self.__coords_yx
            self.__coord = self.__coord_yx

class Selector(object):
    """Makes a selector which stores the selected items in the FeaModel.
    View orientation is also stored.

    Args:
        feamodel (FeaModel): the parent FeaModel object

    Attributes:
        __fea (FeaModel): the parent FeaModel
        __all (bool): If everything is selected all=True
        parts (set): currently selected parts
        areas (set): currently selected areas
        lines (set): currently selected signed lines or arcs
        points (set): currently selected points
        elements (set): currently selected elements
        faces (set): currently selected faces
        nodes (set): currently selected nodes
        orientation (str): view orientation

            'xy': default, x is horizontal, y is vertical
            'yx': y is horizontal, x is vertical

        __parts (set): internal storage set of parts
        __areas (set): internal storage set of areas
        __slines (set): internal storage set of signed lines and arcs
        __points (set): internal storage set of points
        __elements (set): internal storage set of elements
        __faces (set): internal storage set of faces
        __nodes (set): internal storage set of nodes
    """
    def __init__(self, feamodel):
        self.__fea = feamodel
        self.__all = True
        self.__parts = set()
        self.__areas = set()
        self.__slines = set()
        self.__points = set()
        self.__elements = set()
        self.__faces = set()
        self.__nodes = set()
        self.orientation = 'xy'

    @property
    def parts(self):
        """Returns selected parts."""
        if self.__all:
            return self.__fea.parts
        else:
            return self.__parts

    @property
    def areas(self):
        """Returns selected areas."""
        if self.__all:
            return self.__fea.areas
        else:
            return self.__areas

    @property
    def lines(self):
        """Returns selected signed lines or arcs."""
        if self.__all:
            return self.__fea.signlines
        else:
            return self.__slines

    @property
    def points(self):
        """Returns selected points."""
        if self.__all:
            points = [pnt for pnt in self.__fea.points if pnt.arc_center == False]
            return points
        else:
            return self.__points

    @property
    def elements(self):
        """Returns selected elements."""
        if self.__all:
            return self.__fea.elements
        else:
            return self.__elements

    @property
    def faces(self):
        """Returns selected faces."""
        if self.__all:
            return self.__fea.faces
        else:
            return self.__faces

    @property
    def nodes(self):
        """Returns selected nodes."""
        if self.__all:
            return self.__fea.nodes
        else:
            return self.__nodes

    @staticmethod
    def __add_items(inclusive, items, item, parent_list, this_set):
        """Adds an item to items if appropriate.

        Helper function for select_above_all
        """
        if inclusive == False:
            items.add(item)
        else:
            parent_set = set(parent_list)
            if parent_set.issubset(this_set):
                items.add(item)

    def __add_select(self, items):
        """Adds items to the selection set."""
        for item in items:
            if isinstance(item, partmodule.Part):
                self.__parts.add(item)
            elif isinstance(item, geometry.Area):
                self.__areas.add(item)
            elif (isinstance(item, geometry.SignLine)
                  or isinstance(item, geometry.SignArc)):
                self.__slines.add(item)
            elif isinstance(item, geometry.Point):
                self.__points.add(item)
            elif isinstance(item, mesh.Element):
                self.__elements.add(item)
            elif isinstance(item, mesh.Face):
                self.__faces.add(item)
            elif isinstance(item, mesh.Node):
                self.__nodes.add(item)

    def allsel_under(self, sel_type, full=True, byfaces=False):
        """Selects all progeny under current sel_type items.

        Args:
            sel_type (str): the item type to select under

                - 'parts'
                - 'areas'
                - 'lines'
                - 'points' (needs full=True)
                - 'elements'
                - 'faces'

            full (bool): if False selects under one branch

                - solid modeling branch: area->line->keypoints->points
                - mesh branch: elements->faces->nodes

                Cross branch connections:
                - faces are under lines
                - elements are under areas
                - nodes are under lines and points

                If True, selects all items connected to sub items:
                - Parts:

                    - Areas
                    - Elements, element faces + element nodes
                    - Lines, line points

                - Areas:

                    - Elements, element faces + element nodes
                    - Lines, line points

                -Lines:

                    - Faces, nodes
                    - byfaces=True: elements above faces
                    - byfaces=False: elements above line nodes selected

                -Points:

                    - Nodes
                    - byfaces=True: no elements selected
                    - byfaces=False: elements above point nodes selected

            byfaces (bool): If true, the child elements will only be selected
                if a parent face was selected. If false just an element node
                needs to be selected to select the element.

                -This is only applicable when sel_type == 'lines'
        """
        sets = ['parts', 'areas', 'lines', 'points', 'elements', 'faces']
        if isinstance(sel_type, str):
            if sel_type in sets:
                this_set = getattr(self, sel_type)
                for item in this_set:
                    if isinstance(item, partmodule.Part):
                        # areas
                        self.select_below_all('parts')
                        # slines
                        self.select_below_all('areas')
                        # points
                        self.select_below_all('lines')
                        # select related elements if full is on
                        if full:
                            elements = item.elements
                            self.select(elements, also=True)
                            # sel faces
                            self.select_below_all('elements')
                            # sel nodes
                            self.select_below_all('faces')

                    elif isinstance(item, geometry.Area):
                        # slines
                        self.select_below_all('areas')
                        # points
                        self.select_below_all('lines')
                        # select related elements if full is on
                        if full:
                            elements = item.elements
                            self.select(elements, also=True)
                            # sel faces
                            self.select_below_all('elements')
                            # sel nodes
                            self.select_below_all('faces')
                    elif (isinstance(item, geometry.SignLine)
                          or isinstance(item, geometry.SignArc)):
                        # sel points
                        self.select_below_all('lines')
                        # select related faces, nodes, and elements
                        if full:
                            # select faces
                            faces = item.faces
                            self.select(faces, also=True)
                            # sel nodes
                            self.select_below_all('faces')
                            # byfaces logic
                            if byfaces:
                                # sel elements above faces
                                self.select_above_all('faces', inclusive=False)
                            else:
                                # select elements connected to nodes
                                elements = set()
                                for node in self.nodes:
                                    elements.update(node.elements)
                                self.select(elements, also=True)
                    elif isinstance(item, geometry.Point):
                        if full:
                            # sel nodes and elements
                            self.select(item.nodes, also=True)
                            if byfaces == False:
                                for node in item.nodes:
                                    self.select(node.elements, also=True)
            else:
                print('The passed sel_type must be a valid type.')
        else:
            print('The passed sel_type must be a string!')

    def select_below(self):
        """Selects children below currently selected items.

        Only one set must be selected.
        """
        sets = ['parts', 'areas', 'lines', 'points', 'elements', 'faces',
                'nodes']
        sel_sets = []
        for set_name in sets:
            this_set = getattr(self, set_name)
            if len(this_set) > 0:
                sel_sets.append(set_name)
        # check if number of visible sets is one
        if len(sel_sets) == 0 or len(sel_sets) > 1:
            print('Only one set must be selected to use this function.')
        elif len(sel_sets) == 1:
            sel_type = sel_sets[0]
            if sel_type == 'nodes':
                print("One can't select anything below nodes!'")
            else:
                self.select_below_all(sel_type)

    def select_below_all(self, sel_type):
        """Selects children below all currently selected items of sel_type.

        Args:
            sel_type (str): type to select items below

                - 'parts': selects areas
                - 'areas': selects lines (type SignLine and SignArc)
                - 'lines': selects points
                - 'elements': selects faces
                - 'faces': selects nodes
        """
        items = set()
        if sel_type == 'parts':
            for part in self.parts:
                for area in part.areas:
                    items.add(area)
        elif sel_type == 'areas':
            for area in self.areas:
                items.update(area.signlines)
        elif sel_type == 'lines':
            for sline in self.lines:
                items.update(sline.line.points)
        elif sel_type == 'elements':
            for element in self.elements:
                items.update(element.faces)
        elif sel_type == 'faces':
            for face in self.faces:
                items.update(face.nodes)
        self.__add_select(items)

    def select_above(self, inclusive=False):
        """Selects parents above currently selected items.

        Only one set must be selected.
        """
        sets = ['parts', 'areas', 'lines', 'points', 'elements', 'faces',
                'nodes']
        sel_sets = []
        for set_name in sets:
            this_set = getattr(self, set_name)
            if len(this_set) > 0:
                sel_sets.append(set_name)
        # check if number of visible sets is one
        if len(sel_sets) == 0 or len(sel_sets) > 1:
            print('Only one set must be selected to use this function.')
        elif len(sel_sets) == 1:
            sel_type = sel_sets[0]
            if sel_type == 'parts':
                print("One can't select anything above parts!'")
            else:
                self.select_above_all(sel_type, inclusive)

    def select_above_all(self, sel_type, inclusive=False):
        """Selects parents above all currently selected items of sel_type.

        Args:
            sel_type (str): type to select items above

                - 'areas': selects parts
                - 'lines': selects areas
                - 'points': selects signlines or signarcs
                - 'faces': selects elements
                - 'nodes': selects faces

            inclusive (bool): if True, all child entities must be selected for
                the parent to be selected. If false all parents connected to
                the child are selected. Exampes:

                - 'areas' + True: all part areas must be selected to select part
                - 'areas' + False: parts with selected areas in them are selected
                - 'lines' + True: all lines must be selected to select area
                - 'lines' + False: areas with selected lines in them are selected
        """
        items = set()
        if sel_type == 'areas':
            this_set = self.areas
            for area in this_set:
                item = area.parent
                parent_list = area.parent.areas
                self.__add_items(inclusive, items, item, parent_list, this_set)
        elif sel_type == 'lines':
            this_set = self.lines
            for sline in this_set:
                item = sline.parent
                parent_list = sline.parent.signlines
                self.__add_items(inclusive, items, item, parent_list, this_set)
        elif sel_type == 'points':
            this_set = set(self.points)
            for point in this_set:
                for parent_line in point.lines:
                    for sline in parent_line.signlines:
                        item = sline
                        parent_list = sline.points
                        self.__add_items(inclusive, items, item, parent_list,
                                         this_set)
        elif sel_type == 'faces':
            this_set = set(self.faces)
            for face in this_set:
                item = face.element
                parent_list = item.faces
                self.__add_items(inclusive, items, item, parent_list, this_set)
        elif sel_type == 'nodes':
            this_set = set(self.nodes)
            for node in this_set:
                for face in node.faces:
                    item = face
                    parent_list = face.nodes
                    self.__add_items(inclusive, items, item, parent_list,
                                     this_set)
        self.__add_select(items)

    def allselect(self):
        """Selects all items."""
        self.select(items='all')

    def select(self, items='', also=False):
        """Selects an item or a list of items.

        Agrs:
            items (str or item or list): item(s) to select
            also (bool): whether to add the items to the current selection set

                - True: add the items to the current selection set
                - False: only select the passed items
        """
        # empty the selection set if also == False
        if items == '':
            print('You must pass in an item or item(s)!')
        else:
            if also == False:
                self.select_none()
            # Turn on 'all' mode if items == None, otherwise, get the items
            if items == 'all':
                self.__all = True
            else:
                items = self.__fea.get_items(items)
                for (ind, item) in enumerate(items):
                    # if a component was passed, get its child items instead
                    if isinstance(item, components.Component):
                        real_items = item.get_children()
                        items[ind] = real_items
                self.__add_select(items)

    def select_all(self, sel_type='all', also=False):
        """Selects all items of sel_type from the full feamodel.

        All items currently selected items are removed from the selection set
        before the all set is selected.

        Args:
            sel_type (str): the items to select

                - 'all'
                - 'parts'
                - 'areas'
                - 'lines'
                - 'points'
                - 'elements'
                - 'faces'
                - 'nodes'
            also (bool): if False, empties the selection set before selecting
                all items, if True adds new items to the current set
        """
        if also == False:
            self.__all = False
            self.select_none()
        # now select all items of type
        if sel_type == 'all':
            self.select()
        elif sel_type == 'parts':
            self.__add_select(self.__fea.parts)
        elif sel_type == 'areas':
            self.__add_select(self.__fea.areas)
        elif sel_type == 'lines':
            self.__add_select(self.__fea.signlines)
        elif sel_type == 'points':
            points = [pnt for pnt in self.__fea.points if pnt.arc_center == False]
            self.__add_select(points)
        elif sel_type == 'elements':
            self.__add_select(self.__fea.elements)
        elif sel_type == 'faces':
            self.__add_select(self.__fea.faces)
        elif sel_type == 'nodes':
            self.__add_select(self.__fea.nodes)

    def set_orientation(self, orientation):
        """Sets the view orientation to 'xy' or 'yx' order is horiz, vert

        Args:
            orientation (str): the view orientation

                'xy': default, x is horizontal, y is vertical
                'yx': y is horizontal, x is vertical

        """
        self.orientation = orientation

    def deselect(self, items):
        """Removes the passed item or items from the selection set.
        """
        self.__all = False
        for item in items:
            if isinstance(item, partmodule.Part):
                self.__parts.discard(item)
            elif isinstance(item, geometry.Area):
                self.__areas.discard(item)
            elif (isinstance(item, geometry.SignLine)
                  or isinstance(item, geometry.SignArc)):
                self.__slines.discard(item)
            elif isinstance(item, geometry.Point):
                self.__points.discard(item)
            elif isinstance(item, mesh.Element):
                self.__elements.discard(item)
            elif isinstance(item, mesh.Face):
                self.__faces.discard(item)
            elif isinstance(item, mesh.Node):
                self.__nodes.discard(item)

    def deselect_all(self, sel_type):
        """Deselects all items of sel_type.

        All items currently selected items are removed from the selection set
        before the all set is selected.

        Args:
            sel_type (str): the items to select

                - 'all'
                - 'parts'
                - 'areas'
                - 'lines'
                - 'points'
                - 'elements'
                - 'faces'
                - 'nodes'
        """
        if self.__all == True:
            self.__parts = self.parts
            self.__areas = self.areas
            self.__slines = self.lines
            self.__points = self.points
            self.__elements = self.elements
            self.__faces = self.faces
            self.__nodes = self.nodes
            self.__all = False
        if sel_type == 'all':
            self.select_none()
        elif sel_type == 'parts':
            self.__parts = set()
        elif sel_type == 'areas':
            self.__areas = set()
        elif sel_type == 'lines':
            self.__slines = set()
        elif sel_type == 'points':
            self.__points = set()
        elif sel_type == 'elements':
            self.__elements = set()
        elif sel_type == 'faces':
            self.__faces = set()
        elif sel_type == 'nodes':
            self.__nodes = set()

    def select_none(self):
        """Empties out the selection object.
        """
        self.__all = False
        self.__parts = set()
        self.__areas = set()
        self.__slines = set()
        self.__points = set()
        self.__elements = set()
        self.__faces = set()
        self.__nodes = set()

    def print_summary(self):
        """Prints a summary of the number so items selected in each sel_type.
        """
        sets = ['parts', 'areas', 'lines', 'points', 'elements', 'faces',
                'nodes']
        spacer = '----------------------------------'
        print(spacer)
        print("View's currently selected items:")
        for sel_type in sets:
            items = getattr(self, sel_type)
            print(' %s: %i selected' % (sel_type, len(items)))
        print(spacer)

