"""This module stores the Selector object

Selector is used in FeaModel to store the model's selected set.
"""

from . import part
from . import geometry
from . import mesh
from . import components

class Selector(object):
    """Makes a selector which stores the selected items in the FeaModel.

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
            if isinstance(item, part.Part):
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
                    if isinstance(item, part.Part):
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
            for apart in self.parts:
                for area in apart.areas:
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

    def deselect(self, items):
        """Removes the passed item or items from the selection set.
        """
        self.__all = False
        for item in items:
            if isinstance(item, part.Part):
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

