""" This module stores the FeaModel class.

This class is the highest level object in a pycalculix program.
It stores all parts, loads, constraints, mesh, problem, and results_file
objects.
"""

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection  # element plotting
import matplotlib.colors as colors
import matplotlib.cm as cmx
import subprocess # used to launch meshers cgx and gmsh
import os # used to delete files written by cgx
from numpy.core.function_base import linspace # need to make contours

from . import environment
from . import base_classes
from . import geometry
from . import material
from . import components
from . import loads
from . import mesh
from . import part
from . import problem
from . import selector

# element colors, 0-1, 0=black, 1=whate
ECOLOR = '.4'
FCOLOR = '.9'
CMAP = 'jet'

class FeaModel(object):
    """Makes a FeaModel instance.

    Parts, area, lines, arcs, points, nodes, elements, faces, models,
    components, and loads are stored in this object.

    Args:
        model_name (str): project name for this FeaModel, this is a file prefix
        ccx (None or str): path to Calculix ccx solver, pass this when you want
            to overwite the default program location.
            None means the default envionment.CCX is used.
        cgx (None or str): path to Calculix cgx mesher, pass this when you want
            to overwite the default program location.
            None means the default envionment.CGX is used.
        gmsh (None or str): path to gmsh mesher, pass this when you want
            to overwite the default program location.
            None means the default envionment.GMSH is used.

    Attributes:
        fname (str): FeaModel project file name prefix
        points (Itemlist): list of all geometry points
        lines (Itemlist): list of all geometry lines
        areas (Itemlist):list of all geometry areas
        parts (Itemlist): list of all geometry parts
        matls (Itemlist): list of all materials
        components (Itemlist): list of all components
        loads (dict): a dictionary of loads

            | Key (float): the load time point
            | loads[time] = list of loads for that time step
            | See method set_time to change the current time.
            | Time = 0.0 stores constant loads, such as:
            |    material, thickness

        problems (Itemlist): list of problems
        nodes (Meshlist): list of all mesh nodes
        eshape (str): element shape

            - 'quad': quadrilateral elements (Default)
            - 'tri': triangle elements
        eorder (int): element order, 1 or 2

            - 1: elements meshed with corner nodes only
            - 2: (Default) elements meshed with corner and midside nodes
        elements (Meshlist): list of all mesh elements
        faces (list): list of all element faces, includes non-exterior ones
        view (Selector): currently selected items

            Attributes:
                parts: list of selected parts
                areas: list of selected areas
                lines: list of selected signed lines or arcs
                points: list of selected points
                elements: list of selected elements
                faces: list of selected faces
                nodes: list of selected nodes
        time (float): current model time value, defaults to 1.0
        units (dict): a dict to store units for different fields

            Keys:
                - 'displ': displacement or location
                - 'force': force
                - 'stress': stress
                - 'temp': temperature
                - 'density': density (mass/volume)
                - 'time': time

            Returns:
                Text string describing the units for the given key field.

                See the set_units method.

                For example when units have been set to metric, the below
                values are returned.

                - 'dist' --> 'm'
                - 'density' --> 'kg/(m^3)'
                - 'time' --> 's'
    """
    def __init__(self, model_name, ccx=None, cgx=None, gmsh=None):
        self.fname = model_name
        self.points = base_classes.Itemlist()
        self.lines = base_classes.Itemlist()
        self.signlines = base_classes.Itemlist()
        self.areas = base_classes.Itemlist()
        self.parts = base_classes.Itemlist()
        self.matls = base_classes.Itemlist()
        self.components = base_classes.Itemlist()
        self.loads = {}    # store loads by time
        self.problems = base_classes.Itemlist()
        self.nodes = base_classes.Meshlist()
        self.elements = base_classes.Meshlist()
        self.faces = []
        self.view = selector.Selector(self)
        self.time = 1.0 # 0 is model set-up, 1 is first step, etc
        self.eshape = 'quad'
        self.eorder = 2
        self.units = {}

        # fix the paths to the needed programs, ccx, cgx, and gmsh
        if ccx != None:
            environment.CCX = ccx
        if cgx != None:
            environment.CGX = cgx
        if gmsh != None:
            environment.GMSH = gmsh

    def set_ediv(self, items, ediv):
        """Sets the number of elements on the passed line.

        Args:
            items (str or SignLine or SignArc or list): lines to set ediv on

                - str: 'L0'
                - list of str ['L0', 'L1']
                - list of SignLine or SignArc part.bottom or part.hole[-1]

            ediv (int): number of elements to mesh on the line
        """
        items = self.get_items(items)
        for line in items:
            line.set_ediv(ediv)

    def set_units(self, dist_unit='m', cfswitch=False):
        """Sets the units that will be displayed when plotting.

        Picks a unit set based on the passed distance unit.
        That unit set is printed to the console when set.
        Defaults to MKS units (meter-kilogram-second)

        ========  =====  ======  ===========  ===============  ====
        Distance  Force  Stress  Temperature  Density          Time
        ========  =====  ======  ===========  ===============  ====
        'm'       'N'    'Pa'    'K'          'kg/(m^3)'       's'
        'mm'      'N'    'MPa'   'K'          'tonne/(mm^3)'   's'
        'in'      'lbf'  'psi'   'R'          'slinch/(in^3)'  's'
        'ft'      'lbf'  'psf'   'R'          'slug/(ft^3)'    's'
        ========  =====  ======  ===========  ===============  ====

        See get_units method or returning text strings based on unit types.

        Args:
            dist_unit (str): string of distance unit. Options:

                - 'm': meter
                - 'mm': milimeter
                - 'in': inch
                - 'ft': foot

            cfswitch (bool): Celsius/Fahrenheit temperature switch.
                Default is False.

                If True, this switches temperature from K-->C or from R-->F
                Default keeps units as shown in above table.
        """
        keys = ['dist', 'force', 'stress', 'temp', 'density', 'time']

        m_newton = ['m', 'N', 'Pa', 'K', 'kg/(m^3)', 's']
        mm_newton = ['mm', 'N', 'MPa', 'K', 'tonne/(mm^3)', 's']
        in_lbf = ['in', 'lbf', 'psi', 'R', 'slinch/(in^3)', 's']
        ft_lbf = ['ft', 'lbf', 'psf', 'R', 'slug/(ft^3)', 's']
        unit_systems = [m_newton, mm_newton, in_lbf, ft_lbf]
        vals = [usys for usys in unit_systems if usys[0] == dist_unit][0]

        if cfswitch:
            # switches K-->C and R-->F
            newunit = {'K':'C', 'R':'F'}
            vals['temp'] = newunit[vals['temp']]

        # set values
        adict = dict(zip(keys, vals))
        adict['displ'] = adict['dist']
        print('================================================')
        print('Units have been set to %s_%s' % (adict['dist'], adict['force']))
        for key in adict:
            print('For %s use %s' % (key, adict[key]))
        print('================================================')
        self.units = adict

    def get_units(self, *args):
        """Returns units for the passed arguments. Accepts and returns a list
        of strings.

        Options for inputs:
            - 'dist'
            - 'displ' (same units as 'dist')
            - 'force'
            - 'stress'
            - 'temp'
            - 'density'
            - 'time'
        """
        res = []
        for arg in args:
            mystr = ''
            if arg in self.units:
                # check if the requested item is in the units dict
                mystr = ' ('+ self.units[arg] + ')'
            elif arg in base_classes.FIELDTYPE:
                # return units on results dict item, ex: stress, strain etc
                ftype = base_classes.FIELDTYPE[arg]
                if ftype in self.units:
                    mystr = ' ('+ self.units[ftype] + ')'
            res.append(mystr)
        return res

    def set_time(self, time):
        """Sets the time in the FeaModel (preprocessor).

        This time is used when setting loads.

        Args:
            time (float): the time to set
        """
        self.time = time

    def get_item(self, item):
        """Returns an item given a string identifying the item.

        Args:
            item (str): 'A0' or 'P0' or 'L0' etc.
        """
        if item[0] == 'P':
            # get point
            items = self.points
            num = int(item[1:])
            res = [a for a in items if a.id == num]
            return res[0]
        elif item[0] == 'L' or item[1] == 'L':
            # get line
            items = self.signlines
            num = int(item[1:])
            items = [a for a in items if a.get_name() == item]
            return items[0]
        elif item[0] == 'A':
            # get area
            items = self.areas
            num = int(item[1:])
            items = [a for a in items if a.id == num]
            return items[0]
        elif item[0] == 'E':
            # get element
            items = self.elements
            items = [a for a in items if a.get_name() == item]
            return items[0]
        elif item[0] == 'N':
            # get node
            items = self.nodes
            items = [a for a in items if a.get_name() == item]
            return items[0]
        else:
            print('Unknown item! Please pass the name of a point, line or area!')

    def get_items(self, items):
        """Returns a list of correctly typed items.

        Input can be a single string or item, or a list of strings identifying
        items.
        """
        res_items = base_classes.listify(items)
        for ind, item in enumerate(res_items):
            if isinstance(item, str):
                res_items[ind] = self.get_item(item)
        return res_items

    def make_matl(self, name):
        """Makes and returns a new material.

        Args:
            name (str): the material's name
        """
        mat = material.Material(name)
        self.matls.append(mat)
        return mat

    def make_part(self):
        """Makes and returns a new part."""
        #p = part.Part(self)
        #self.parts.append(p)
        return part.Part(self)

    def make_problem(self, problem_type='struct', parts='all'):
        """Makes and returns a new problem, which can be solved.

        Args:
            problem_type (str): problem type, options:
                'struct': structural
            parts (str Part or list of Part): Parts the model will analyze.

                Options:
                    - 'all': add all parts. This is the default.
                    - Part: add the single part
                    - list of Part: add these parts
        """
        if parts == 'all':
            parts = self.parts
        prob = problem.Problem(self, problem_type, parts)
        return prob

    def plot_elements(self, fname='', display=True, title='Elements',
                      enum=False, nshow=False, nnum=False):
        """Plots the selected elements.

        Args:
            fname (str): png image file prefix, if given, image will be saved
            display (bool): if True, interactive plot will be shown, if False
                plot will not be shown
            title (str): the plot's title
            enum (bool): if True element numbers are plotted
            nshow (bool): True=plot nodes, False=don't plot them
            nnum (bool): if True node numbers are plotted
        """
        nodes = self.view.nodes
        elements = self.view.elements
        if len(elements) > 0:
            # plotting elements
            fig, ax = plt.subplots()
            polys = []
            for element in elements:
                poly = element.get_poly()
                polys.append(poly)
            coll = PatchCollection(polys, facecolors=FCOLOR, edgecolors=ECOLOR)
            ax.add_collection(coll)

            # plot element numbers
            if enum:
                for element in elements:
                    element.label(ax)

            # plot nodes, this is quicker than individual plotting
            if nshow:
                axs = [node.y for node in nodes]
                rads = [node.x for node in nodes]
                ax.scatter(axs, rads, s=7, color='black')
                if nnum:
                    for node in nodes:
                        node.label(ax)

            # set units
            [d_unit] = self.get_units('dist')
            plt.title(title)
            plt.xlabel('axial, y'+d_unit)
            plt.ylabel('radial, x'+d_unit)
            plt.axis('scaled')

            # extract nodes for elements if no nodes are selected
            if len(nodes) == 0:
                tmp = set()
                for element in elements:
                    tmp.update(element.nodes)
                nodes = list(tmp)
            # extract max and min for plot window
            radials = [n.x for n in nodes]
            axials = [n.y for n in nodes]

            # finish pot
            base_classes.plot_set_bounds(plt, axials, radials)
            base_classes.plot_finish(plt, fname, display)

        else:
            # no elements exist or no elemnts are selected
            res = ''
            if len(self.elements) == 0:
                res = 'No elements exist! Try meshing your parts!'
            else:
                res = 'No elements are selected! Select some!'
            print(res)

    def plot_pressures(self, fname='', display=True):
        """Plots the load step pressures on the passed part(s) or selected
        elements.

        This is an element plot, with arrows showing pressure magnitude and
        directions.

        Args:
            fname (str): png image file prefix, if given, image will be saved
            display (bool): if True, interactive plot will be shown, if False
                plot will not be shown. Default = True
        """
        elements = self.view.elements
        faces = self.view.faces
        nodes = self.view.nodes

        # plot all elements and store length, length determins min pressure arrow
        if len(elements) > 0:
            # plotting elements
            fig = plt.figure()
            ax = fig.add_subplot(111, aspect='equal')

            # plot polys and calculate avg face length
            polys = []
            face_len = []
            for element in elements:
                face_len.append(element.face[1].length())
                poly = element.get_poly()
                polys.append(poly)
            coll = PatchCollection(polys, edgecolors=ECOLOR, facecolors=FCOLOR)
            ax.add_collection(coll)
            face_len = sum(face_len)/len(face_len)

            # store pressures we'll want to plot: list of [face, pval]
            plist = []
            for load in self.loads[self.time]:
                if load.ltype in  ['press', 'press_fluid']:
                    loadlist = load.get_list()
                    for [face, pval] in loadlist:
                        if face in faces:
                            plist.append([face, pval])
            pressures = [abs(pval) for [face, pval] in plist]

            # set max and min pressure bounds
            pmin = min(pressures)
            pmax = max(pressures)
            # extract nodes for elements if no nodes are selected
            if len(nodes) == 0:
                tmp = set()
                for element in elements:
                    tmp.update(element.nodes)
                nodes = list(tmp)
            radials = [p.x for p in nodes]
            axials = [p.y for p in nodes]
            # arrow length = arrow_min + (pval-pmin)*mult
            arrow_min = face_len
            mult = 0
            if pmax != pmin:
                adelta = max(axials) - min(axials)
                rdelta = max(radials) - min(radials)
                delta = max(adelta, rdelta)
                dist = delta*0.2
                mult = dist/(pmax - pmin)

            # make tick list for later plot, and color map
            tick_list = [pmin]  # default to plot one val
            cmap = colors.ListedColormap(['b', 'b']) # default to plot one val
            if pmax != pmin:
                # we have a range of values we're plotting
                tick_list = linspace(pmin, pmax, 8)
                cmap = plt.get_cmap(CMAP)

            # set color contours for arrows
            cnorm = colors.Normalize(vmin=pmin, vmax=pmax)
            scalarmap = cmx.ScalarMappable(norm=cnorm, cmap=cmap)
            scalarmap.set_array([])

            # make arrows, store axials + radials
            for [face, pval] in plist:
                [face_point, unit] = face.get_mnorm()
                arrow_length = arrow_min + (abs(pval) - pmin)*mult
                other_point = face_point + unit*arrow_length
                radials.append(other_point.x)
                axials.append(other_point.y)

                # compression
                tail = other_point
                delta = face_point - other_point
                if pval < 0:
                    # tension
                    tail = face_point
                    delta = other_point - face_point

                hw = face_len*0.2
                hl = face_len*0.3
                colorVal = scalarmap.to_rgba(abs(pval))
                plt.arrow(tail.y, tail.x, delta.y, delta.x,
                          color=colorVal, head_width=hw, head_length=hl,
                          length_includes_head=True)

            # set the horizontal and vertical axes
            base_classes.plot_set_bounds(plt, axials, radials)

            # set units and titles
            [d_unit, p_unit, t_unit] = self.get_units('dist', 'stress', 'time')
            tstr = 'Pressures %s\nTime=%f%s' % (p_unit, self.time, t_unit)
            plt.title(tstr)
            plt.xlabel('axial, y'+d_unit)
            plt.ylabel('radial, x'+d_unit)

            # set the colorbar
            plt.colorbar(scalarmap, orientation='vertical', ticks=tick_list)
            base_classes.plot_finish(plt, fname, display)

        else:
            res = ''
            if len(self.elements) == 0:
                res = 'No elements exist! Try meshing your parts!'
            else:
                res = 'No elements are selected! Select some!'
            print(res)

    def plot_constraints(self, fname='', display=True):
        """Plots the constraints on the passed part(s) or selected items.

        This is an element and node plot, with arrows showing displacement
        constraints.

        Args:
            fname (str): png image file prefix, if given, image will be saved
            display (bool): if True, interactive plot will be shown, if False
                plot will not be shown. Default = True
        """
        elements = self.view.elements
        nodes = self.view.nodes

        # plot all elements
        if len(elements) > 0:
            # plotting elements
            fig = plt.figure()
            ax = fig.add_subplot(111, aspect='equal')

            # plot polys and calculate avg face length
            polys = []
            face_len = []
            for element in elements:
                face_len.append(element.face[1].length())
                poly = element.get_poly()
                polys.append(poly)
            coll = PatchCollection(polys, edgecolors=ECOLOR, facecolors=FCOLOR)
            ax.add_collection(coll)
            face_len = sum(face_len)/len(face_len)

            # store displacements we'll plot: list of [node, dict ux,uy,uz]
            ulist = []
            vals = []
            for load in self.loads[self.time]:
                if load.ltype in  ['ux', 'uy', 'uz']:
                    alist = load.get_list()
                    for nodelist in alist:
                        node = nodelist[0]
                        if node in nodes:
                            ulist += [nodelist]
                            vals += list(nodelist[1].values())

            # check min and max bounds
            pmin = min(vals)
            pmax = max(vals)

            # make tick list for later plot, and color map
            tick_list = [pmin]  # default to plot one val
            cmap = colors.ListedColormap(['b', 'b']) # default to plot one val
            if pmax != pmin:
                # we have a range of values we're plotting
                tick_list = linspace(pmin, pmax, 8)
                cmap = plt.get_cmap(CMAP)

            # set color contours for arrows
            cnorm = colors.Normalize(vmin=pmin, vmax=pmax)
            scalarmap = cmx.ScalarMappable(norm=cnorm, cmap=cmap)
            scalarmap.set_array([])

            # make arrows for displacements
            # extract nodes for elements if no nodes are selected
            if len(nodes) == 0:
                tmp = set()
                for element in elements:
                    tmp.update(element.nodes)
                nodes = list(tmp)
            radials = [p.x for p in nodes]
            axials = [p.y for p in nodes]
            pvect = {'ux':geometry.Point(1, 0, 0),
                     'uy':geometry.Point(0, 1, 0),
                     'uz':geometry.Point(0, 0, 1)}
            for [node, udict] in ulist:
                for (key, val) in udict.items():
                    headw = face_len*0.4
                    headl = headw
                    point = geometry.Point(node.x, node.y, node.z)
                    unit = pvect[key]

                    # nonzero displacement, draw from node
                    tail = point
                    head = tail + unit*val
                    store_point = head
                    if val == 0:
                        # zero displ, draw to node
                        head = point
                        tail = head - unit*face_len
                        store_point = tail
                    delta = head - tail
                    radials.append(store_point.x)
                    axials.append(store_point.y)

                    colorVal = scalarmap.to_rgba(val)
                    plt.arrow(tail.y, tail.x, delta.y, delta.x,
                              color=colorVal, head_width=headw,
                              head_length=headl, length_includes_head=True)

            # set the horizontal and vertical axes
            base_classes.plot_set_bounds(plt, axials, radials)

            # set units + titles
            [d_unit, t_unit] = self.get_units('dist', 'time')
            tstr = 'Constraints %s\nTime=%f%s' % (d_unit, self.time, t_unit)
            plt.title(tstr)
            plt.xlabel('axial, y'+d_unit)
            plt.ylabel('radial, x'+d_unit)

            # set the colorbar
            plt.colorbar(scalarmap, orientation='vertical', ticks=tick_list)
            base_classes.plot_finish(plt, fname, display)

        else:
            res = ''
            if len(self.elements) == 0:
                res = 'No elements exist! Try meshing your parts!'
            else:
                res = 'No elements are selected! Select some!'
            print(res)

    def plot_geometry(self, fname='', display=True, title='Geometry', pnum=True,
                      lnum=True, anum=True, afill=True):
        """Plots selected geometry items: areas, lines, points.

        Defaults to displaying labels on: areas, lines, and points, and filling
        in the areas.

        Args:
            fname (str): png image file prefix, if given, image will be saved
            display (bool): if True, interactive plot will be shown, if False
                plot will not be shown. Default = True
            title (str): the plot's title
            pnum (bool): if True all Point objects will be labeled
            lnum (bool): if True all Line and Arc objects will be labeled
            anum (bool): if True all Area objects will be labeled
            afill (bool): if True all Area objects will be filled
        """
        # this method plots the part
        # check out: http://nickcharlton.net/posts/drawing-animating-shapes-matplotlib.html

        points = self.view.points
        lines = self.view.lines
        areas = self.view.areas

        # start plotting
        fig = plt.figure()
        ax = fig.add_subplot(111)

        #plot areas
        for area in areas:
            area.plot(ax, anum, afill)

        #plot lines
        for sline in lines:
            sline.plot(ax, lnum)

        # plot points
        for point in points:
            point.plot(ax, pnum)

        # set the horizontal and vertical axes
        radials = [point.x for point in points]
        axials = [point.y for point in points]
        base_classes.plot_set_bounds(plt, axials, radials)

        # set units
        [d_unit] = self.get_units('dist')

        # show plot
        plt.title(title)
        plt.xlabel('axial, y'+d_unit)
        plt.ylabel('radial, x'+d_unit)
        ax.set_aspect('equal')
        base_classes.plot_finish(plt, fname, display)

    def __get_cname(self, items):
        """Returns a component name prefix, for labeling lists of items.

        Args:
            items (list): a list of the items we'll be making a component of
        """
        cname = ''
        if len(items) == 1:
            if items[0] == 'all':
                cname = 'all'
            else:
                cname = items[0].get_name()
        else:
            cname = items[0].get_name()+'-'+items[-1].get_name()
        return cname

    def __get_make_comp(self, comp):
        """Stores component if it doesn't exist, returns it if it does.

        Args:
            comp (Component): component we want to store/get
        """
        if comp not in self.components:
            comp = self.components.append(comp)
        else:
            ind = self.components.index(comp)
            comp = self.components[ind]
        return comp

    def __add_load(self, load, time):
        """Adds a load to the FeaModel.

        Args:
            load (Load or Load_Linear): th eload to add
            time (float): the load's time point. The first time point is 1.0
        """
        if time in self.loads:
            self.loads[time].append(load)
        else:
            self.loads[time] = [load]
        return load

    def set_gravity(self, grav, items):
        """Sets gravity on the elements in items.

        Assumes gravity acts in the -x direction with magnitude grav.

        Args:
            grav (float): gravity acceleration, MUST BE POSTIVE
            items (Area or Part or list): items gravity acts on

                - Area: gravity acts on elements in this area
                - Part: gravity acts on elements in this part
                - list of Part or Area: gravity acts on their elements
        """
        items = self.get_items(items)
        ctype = 'elements'
        cname = self.__get_cname(items)

        # make compoenet
        comp = components.Component(items, ctype, cname)
        comp = self.__get_make_comp(comp)

        # make load
        ltype = 'gravity'
        time = self.time
        load = loads.ConstLoad(ltype, comp, grav)
        self.__add_load(load, time)
        return load

    def set_rpm(self, rpm, items):
        """Sets rpm rotation load on the items.

        Args:
            rpm (float): rotation per minute
            items (Area or Part or list): items rotation acts on

                - Area: rotation acts on elements in this area
                - Part: rotation acts on elements in this part
                - list of Part or Area: rotation acts on their elements
        """
        # applies rpm to the items
        items = self.get_items(items)
        ctype = 'elements'
        cname = self.__get_cname(items)

        # make compoenet
        comp = components.Component(items, ctype, cname)
        comp = self.__get_make_comp(comp)

        # make load
        ltype = 'rpm'
        time = self.time
        load = loads.ConstLoad(ltype, comp, rpm)
        self.__add_load(load, time)
        return load

    def set_radps(self, radps, items):
        """Sets radians per second rotation load on the items.

        Args:
            radps (float): radians per second
            items (Area or Part or list): items rotation acts on

                - Area: rotation acts on elements in this area
                - Part: rotation acts on elements in this part
                - list of Part or Area: rotation acts on their elements
        """
        items = self.get_items(items)
        ctype = 'elements'
        cname = self.__get_cname(items)

        # make compoenet
        comp = components.Component(items, ctype, cname)
        comp = self.__get_make_comp(comp)

        # make load
        ltype = 'radps'
        time = self.time
        load = loads.ConstLoad(ltype, comp, radps)
        self.__add_load(load, time)
        return load

    def set_fluid_press(self, items, rho, g, xo, po):
        """This sets a fluid presure load on items.

        This fluid pressure is dependednt on the x axis.
        g must be positive.

            - P = f(x)
            - P(xo) = po
            - P(x) = po + rho*g*(xo - x)

        Args:
            items (str or Line or Arc or list): items to set pressure on

                - str: string name of Line or Arc item, for example 'L0'
                - Line or Arc: set presssure on this
                - list or Line or Arc: set pressure on these
            rho (float): fluid density in mass/volume
            g (+float): gravity in dist/(t^2), MUST BE POSITIVE
            xo (float): sea level height, MUST BE POSITIVE
            po (float): sea level pressure, MUST BE POSITIVE
        """
        items = self.get_items(items)
        ctype = 'faces'
        cname = self.__get_cname(items)

        # make compoenet
        comp = components.Component(items, ctype, cname)
        comp = self.__get_make_comp(comp)

        # make load
        ltype = 'press_fluid'
        mult = rho*g
        load = loads.LinearLoad(ltype, comp, po, mult, xo)
        self.__add_load(load, self.time)
        return load

    def set_load(self, ltype, items, lval, ldir=None):
        """Sets a pressure or force load on line(s).

        The force load is divided by the number of nodes and applied on
        each node.

        The pressure load is applied to the child faces under the line.
        Positive is compression, negative is tension.

        Args:
            ltype (str): 'press' or 'force'
            items (str or Line or Arc or list): items to set load on

                - str: string name of Line or Arc item, for example 'L0'
                - Line or Arc: set load on this
                - list or Line or Arc: set load on these
            lval (float): load value.

                - For ltype = 'press' this is in force/area units
                - For ltype = 'force' this is in force units
            ldir (None or str): load direction. Defaults to None

                - str: when ltype='load', we need to set ldir to 'x' or 'y'
                - None: when ltype='press'
        """
        items = self.get_items(items)
        ctype = 'nodes'
        if ltype == 'press':
            ctype = 'faces'
        cname = self.__get_cname(items)

        # make compoenet
        comp = components.Component(items, ctype, cname)
        comp = self.__get_make_comp(comp)

        # make load
        if ltype == 'force':
            ltype = 'f'+ldir # for example fx
        load = loads.ConstLoad(ltype, comp, lval)
        self.__add_load(load, self.time)
        return load

    def set_constr(self, ltype, items, axis, val=0.0):
        """Sets a displacement constraint on the passed item(s).

        Args:
            ltype (str): 'fix' or 'displ'

                - 'fix': val arg should not be passed
                - 'displ': val arg must be passed
            items (str, SignLine, SignArc, Point or list): item(s) to apply the
                constraint on

                - str: this string must be a line or point name
                - SignLine or SignArc: the constraint will be applied
                    to this item
                - Point: the constraint will be applied on this item
                - list of SignLine or SignArc: the constraint will be appplied
                    on these
                - list of Point: the constraint will be appplied on these
                - list of str names: pass line or point names, constr
                    applied on them

            axis (str): load axis, 'x' or 'y'
            val (float): displacement value, defaults to 0.0, only needs to be
                set when ltype = 'displ'
        """
        items = self.get_items(items)
        cname = self.__get_cname(items)
        ctype = 'nodes'

        # make compoenet
        comp = components.Component(items, ctype, cname)
        comp = self.__get_make_comp(comp)

        # make load
        ltype = 'u'+axis # for example ux
        load = loads.ConstLoad(ltype, comp, val)
        self.__add_load(load, self.time)
        return load

    def set_eshape(self, eshape='quad', eorder=2):
        """Sets the element shape and order to use when meshing the model.

        Args:
            eshape (str): element shape

                - 'quad': quadrilaterials
                - 'tri': triangles
            eorder (int): element order, default=2

                - 1: corner nodes only (3 and 4 noded elements)
                - 2: corder nodes and midside nodes (6 and 8 noded elements)
        """
        self.eshape = eshape # quad or tri
        self.eorder = eorder # 1 or 2

    def set_etype(self, etype, items, thick=None):
        """Sets the element type, and thickness on areas or parts.

        Args:
            etype (str): element type

                - 'plstress': plane stress
                - 'plstrain': plane strain
                - 'axisym': axisymmetric
            items (str or Area or Part or list): set element type on these

                - str: string name of Area or Part item. Example: 'A0', 'PART0'
                - Area: set element type on the elements in this area
                - Part: set element type on the elements in this part
                - list of Area or Part: set element type on their elements
            thick (float or None): element thickness

                - None: default, used for axisymmetric element type
                - float: thickness value to use for plane stress or plane
                    strain elements
        """
        items = self.get_items(items)

        # set the element types on the areas, this is used to fix
        # elements when importing them from the inp file
        for item in items:
            if isinstance(item, geometry.Area):
                item.set_etype(etype)
            if isinstance(item, part.Part):
                for area in item.areas:
                    area.set_etype(etype)

        # set a thickness component if needed
        if etype != 'axisym' and thick != None:

            # make component for element thickness
            cname = self.__get_cname(items)
            ctype = 'nodes'
            comp = components.Component(items, ctype, cname)
            comp = self.__get_make_comp(comp)

            # add load
            ltype = 'nodal_thickness'
            time = 0.0
            load = loads.ConstLoad(ltype, comp, thick)
            self.__add_load(load, time)
            return load
        return None

    def set_matl(self, matl, items):
        """Sets the material on Part or Area items.

        Args:
            matl (Material): material to assign to items
            items (Part or Area or list): items that will have matl assigned

                - Part: assign matl to this
                - Area: assign matl to this
                - list of Area or Part: assign matl to these
        """
        items = self.get_items(items)
        cname = self.__get_cname(items)

        # store the material if it's not already
        if matl not in self.matls:
            matl = self.matls.append(matl)

        # make component
        ctype = 'elements'
        comp = components.Component(items, ctype, cname)
        comp = self.__get_make_comp(comp)

        # make load
        ltype = 'matl'
        time = 0.0
        load = loads.ConstLoad(ltype, comp, matl)
        self.__add_load(load, time)
        return load

    def __read_inp(self, fname):
        """Reads in the mesh from a Calculix inp file.

        All nodes, elements, and faces are read. Child nodes, elements and
        faces are assigned to their geometry parents (points, lines, arcs,
        areas, parts.

        Args:
            fname (str): file name to read, must include '.inp' extension
        """
        f = open(fname, 'r')
        mode = None
        set_name = None
        set_type = None

        items = [] # holder for nodes or elements in nsets or esets
        N = base_classes.Meshlist() # store nodes
        E = base_classes.Meshlist() # store elements, allows renumbering before putting int model
        F = [] # store faces
        sets = {'E':{}, 'N':{}} # store sets
        etype = ''

        # read in input file
        for line in f:
            if line[0] != '*':
                if mode == 'nmake':
                    L = line.split(',')
                    L = [a.strip() for a in L]
                    (nnum, x, y, z) = (int(L[0]), float(L[1]), float(L[2]), float(L[3]))
                    node = mesh.Node(nnum, x, y, z)
                    N.append(node)
                elif mode == 'emake':
                    L = line.split(',')
                    L = [int(a.strip()) for a in L]
                    enum = L[0]
                    nlist = [N.idget(a) for a in L[1:]]
                    e = mesh.Element(enum, etype, nlist)
                    faces = e.faces
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
                    todel.append({'set_type':set_type, 'set_name':set_name})
        # delete the empty sets
        for adict in todel:
            (set_type, set_name) = (adict['set_type'], adict['set_name'])
            del sets[set_type][set_name]
            #print('Empty set type:%s name:%s deleted' % (set_type, set_name))

        # this resets the min element to number 1
        if E.get_minid() > 1:
            E.set_minid(1)

        #-----------------------------------
        # Node and element assignment back onto parts, areas, lines, points
        #-----------------------------------
        # assign elements + nodes to feamodel
        self.elements = E
        self.faces = F

        # remove arc center ndoes from imported node set
        # those nodes have no elements under them
        torem = []
        for node in N:
            if len(node.elements) == 0:
                torem.append(node)
        for node in torem:
            N.remove(node)
        # assign nodes to feamodel
        self.nodes = N

        for apart in self.parts:
            # assign part element and node sets
            pname = apart.get_name()
            apart.elements = sets['E'][pname]
            apart.nodes = sets['N'][pname]

            # assign all nodes and elements to areas, fix the element types
            for area in apart.areas:
                aname = area.get_name()
                area.elements = sets['E'][aname]
                area.nodes = sets['N'][aname]
                area.set_child_ccxtypes() #paint element types on elements

            # assign the child nodes to points
            pts = apart.points
            for pt in pts:
                ndist = []
                for n in apart.nodes:
                    p_tmp = geometry.Point(n.x, n.y)
                    p_tmp = pt - p_tmp
                    dist = p_tmp.length()
                    ndist.append({'dist':dist, 'node':n})
                # sort the list by dist, sorts low to high
                ndist = sorted(ndist, key=lambda k: k['dist'])
                pt.nodes = [ndist[0]['node']]
                #print('Point %s = node %s' % (pt, pt.nodes))

            # assign the nodes and n1 and faces to lines
            slines = apart.signlines
            for sline in slines:
                lname = sline.get_name()
                area = sline.parent
                nodes = sets['N'][sline.line.get_name()]
                n1 = [n for n in nodes if n.order == 1]
                sline.nodes = nodes
                sline.n1 = n1
                # make a set of all faces that contain the line's node
                allfaces = set()
                for node in nodes:
                    allfaces.update(node.faces)
                # make a set of all faces on the line only
                linefaces = set()
                for face in allfaces:
                    if set(face.nodes).issubset(set(nodes)):
                        linefaces.add(face)
                # and only the ones whose element centers are in the area
                faces = set()
                for face in linefaces:
                    is_lineface = area.contains_point(face.element.center)
                    if is_lineface:
                        faces.add(face)
                sline.faces = list(faces)
        print('Elements: %i' % len(E))
        print('Nodes: %i' % len(N))
        print('Done reading Calculix/Abaqus .inp file')

    def mesh(self, fineness=1.0, mesher='gmsh'):
        """Meshes all parts.

        Args:
            fineness (float): 0.0001 - 1.0, how fine the mesh is.

                - Low numbers are very fine, higher numbers are coarser.
            mesher (str): the mesher to use
            
                - 'gmsh': mesh with Gmsh, this is reccomended, it allows holes
                - 'cgx': mesh with Calculix cgx, it doesn't allow holes
        """
        if mesher == 'gmsh':
            self.__mesh_gmsh(fineness)
        elif mesher == 'cgx':
            self.__mesh_cgx(fineness)

    def __mesh_gmsh(self, fineness):
        """Meshes all parts using the Gmsh mesher.

        Args:
            fineness (float): 0.0001 - 1.0, how fine the mesh is.

                Low numbers are very fine, higher numbers are coarser.
        """
        geo = []
        ids = {}
        ids['line'] = {}
        ids['plane_surface'] = {}

        # write all points
        for pt in self.points:
            linestr = 'Point(%i) = {%f, %f, %f};' % (pt.id, pt.x, pt.y, 0.0)
            geo.append(linestr)

        # start storing an index number
        ind = self.points[-1].id + 1

        # write all lines
        for line in self.lines:
            ln = line.id
            ids['line'][ln] = ind
            p1 = line.pt(0).id
            p2 = line.pt(1).id
            linestr = ''
            if isinstance(line, geometry.Arc):
                # line is arc
                pc = line.actr.id
                linestr = 'Circle(%i) = {%i, %i, %i};' % (ind, p1, pc, p2)
            else:
                # straight line
                linestr = 'Line(%i) = {%i,%i};' % (ind, p1, p2)
            geo.append(linestr)

            # set division if we have it
            if line.ediv != None:
                ndiv = line.ediv+1
                esize = line.length()/line.ediv
                if self.eshape == 'quad':
                    ndiv = line.ediv/2+1
                    esize = esize*2
                    # this is needed because quad recombine
                    # splits 1 element into 2
                linestr = 'Transfinite Line{%i} = %i;' % (ind, ndiv)
                print('LINE ELEMENT SIZE: %f, MAKES %i ELEMENTS'
                      % (line.length()/line.ediv, line.ediv))
                geo.append(linestr)
                geo.append('Characteristic Length {%i,%i} = %f;'
                           % (p1, p2, esize))
            ind += 1

        # write all areas
        for area in self.areas:
            if area.closed:
                aid = area.id
                aname = area.get_name()
                loop_ids = []
                loops = [area.signlines] + area.holes
                for loop in loops:
                    linestr = 'Line Loop(%i) = ' % (ind)
                    loop_ids.append(str(ind))
                    line_ids = []
                    for line in loop:
                        lid = ids['line'][line.line.id]
                        prefix = ''
                        if line.sign == -1:
                            prefix = '-'
                        line_ids.append('%s%i' % (prefix, lid))
                    linestr = linestr + '{'+ ','.join(line_ids)+'};'
                    geo.append(linestr)
                    ind += 1
                loop_ids = ','.join(loop_ids)
                geo.append('Plane Surface(%i) = {%s};' % (ind, loop_ids))
                ids['plane_surface'][aid] = ind
                geo.append("Physical Surface('%s') = {%i};" % (aname, ind))
                ind += 1

        # write part area components
        for p in self.parts:
            # make components for each part
            line = "Physical Surface('%s') = " % (p.get_name())
            area_ids = []
            for area in p.areas:
                if area.closed:
                    aid = ids['plane_surface'][area.id]
                    area_ids.append(str(aid))
            line = line + '{' + ','.join(area_ids) + '};'
            geo.append(line)

        # write all line componenets so we can get nodes out
        for L in self.lines:
            lid = ids['line'][L.id]
            line = "Physical Line('%s') = {%i};" % (L.get_name(), lid)
            geo.append(line)

        # write node componenets
        # node list is not produced by gmsh
        for pt in self.points:
            linestr = "Physical Point('%s') = {%i};" % (pt.get_name(), pt.id)
            geo.append(linestr)

        # set the meshing options
        geo.append('Mesh.CharacteristicLengthFactor = '+str(fineness)+
                   '; //mesh fineness')
        geo.append('Mesh.RecombinationAlgorithm = 1; //blossom')
        #geo.append('Mesh.Lloyd = 1; //smoothing algorithm')

        if self.eshape == 'quad':
            geo.append('Mesh.RecombineAll = 1; //turns on quads')
            geo.append('Mesh.SubdivisionAlgorithm = 1; // quadrangles only')
            #geo.append('Mesh.RecombinationAlgorithm = 1; //turns on blossom needed for quad')

        eo = self.eorder
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
        fname = self.fname+'.geo'
        fout = self.fname+'.inp'
        f = open(fname, 'w')
        for line in geo:
            #print (line)
            f.write(line+'\n')
        f.close()
        print('File: %s was written' % fname)

        # run file in bg mode, -2 is 2d mesh
        astr = "%s %s -2 -o %s" % (environment.GMSH, fname, fout)
        print(astr)
        subprocess.call(astr, shell=True)
        print('File: %s was written' % fout)
        print('Meshing done!')

        # write gmsh msh file
        astr = "%s %s -2 -o %s" % (environment.GMSH, fname, self.fname+'.msh')
        subprocess.call(astr, shell=True)
        print('File: %s.msh was written' % self.fname)

        # read in the calculix mesh
        self.__read_inp(self.fname+'.inp')

    def __mesh_cgx(self, fineness):
        """Meshes all parts using the Calculix cgx mesher.

        Args:
            fineness (float): 0.0001 - 1.0, how fine the mesh is.

                Low numbers are very fine, higher numbers are coarser.
        """
        fbd = []
        comps = []
        cfiles = []

        # Calculix CGX elements
        # axisymmetric
        cgx_elements = {}
        cgx_elements['tri2axisym'] = 'tr6c'
        cgx_elements['tri1axisym'] = 'tr3c'
        cgx_elements['quad2axisym'] = 'qu8c'
        cgx_elements['quad1axisym'] = 'qu4c'
        # plane stress
        cgx_elements['tri2plstress'] = 'tr6s'
        cgx_elements['tri1plstress'] = 'tr3s'
        cgx_elements['quad2plstress'] = 'qu8s'
        cgx_elements['quad1plstress'] = 'qu4s'
        # plane strain
        cgx_elements['tri2plstrain'] = 'tr6e'
        cgx_elements['tri1plstrain'] = 'tr3e'
        cgx_elements['quad2plstrain'] = 'qu8e'
        cgx_elements['quad1plstrain'] = 'qu4e'

        num = 1.0/fineness
        emult = int(round(num)) # this converts fineness to mesh multiplier

        # write all points
        for pt in self.points:
            linestr = 'pnt %s %f %f %f' % (pt.get_name(), pt.x, pt.y, 0.0)
            fbd.append(linestr)
            # gmsh can't make node componenets so don't do it in cgx
            #L = 'seta %s p %s' % (pt.get_name(), pt.get_name())
            #comps.append(L)

        # write all lines
        for line in self.lines:
            lname = line.get_name()
            pt1 = line.pt(0).get_name()
            pt2 = line.pt(1).get_name()
            linestr = ''
            if isinstance(line, geometry.Arc):
                # line is arc
                pctr = line.actr.get_name()
                linestr = 'line %s %s %s %s' % (lname, pt1, pt2, pctr)
            else:
                # straight line
                linestr = 'line %s %s %s' % (lname, pt1, pt2)
            # set division if we have it
            if line.ediv != None:
                ndiv = self.eorder*line.ediv
                linestr += ' '+str(int(ndiv))
            fbd.append(linestr)
            linestr = 'seta %s l %s' % (lname, lname)
            comps.append(linestr)
            cfiles.append(lname)

        # write all areas
        for area in self.areas:
            if area.closed:
                linestr = 'gsur '+area.get_name()+' + BLEND '
                line_ids = []
                for line in area.signlines:
                    lname = '+ '+line.get_name()
                    if line.sign == -1:
                        lname = '- '+line.get_name()[1:]
                    line_ids.append(lname)
                linestr = linestr + ' '.join(line_ids)
                fbd.append(linestr)
                # add area component, nodes + elements
                astr = 'seta %s s %s' % (area.get_name(), area.get_name())
                fbd.append(astr)
                cfiles.append(area.get_name())

        # write part area components
        for apart in self.parts:
            # make components for each part
            # seta P0 s s0 s1
            line = 'seta %s s ' % apart.get_name()
            cfiles.append(apart.get_name())
            area_ids = []
            for area in apart.areas:
                if area.closed:
                    area_ids.append(area.get_name())
            line = line + ' '.join(area_ids)
            fbd.append(line)

            # mesh all areas
            for area in apart.areas:
                aname = area.get_name()
                estr = self.eshape+str(self.eorder)+area.etype
                cgx_etype = cgx_elements[estr]
                fbd.append('elty %s %s' % (aname, cgx_etype))
            fbd.append('div all mult %i' % emult)
            fbd.append('mesh all')

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
        fname = self.fname+'.fbd'
        f = open(fname, 'w')
        for line in fbd:
            #print (line)
            f.write(line+'\n')
        f.close()
        print('File: %s was written' % fname)

        # run file in bg mode
        p = subprocess.call("%s -bg %s" % (environment.CGX, fname), shell=True)
        print('Meshing done!')

        # assemble the output files into a ccx input file
        inp = []
        files = ['all.msh']
        files += [f+'.nam' for f in cfiles]
        for fname in files:
            infile = open(fname, 'r')
            for line in infile:
                # cgx adds E and N prfixes on sets after =, get rid of these
                if '=' in line and fname != 'all.msh':
                    L = line.split('=')
                    line = L[0] + '=' + L[1][1:]
                    inp.append(line.strip())
                else:
                    inp.append(line.strip())
            infile.close()

            # delete file
            os.remove(fname)

        # write out inp file
        fname = self.fname+'.inp'
        outfile = open(fname, 'w')
        for line in inp:
            #print (line)
            outfile.write(line+'\n')
        outfile.close()
        print('File: %s was written' % fname)

        # read in the calculix mesh
        self.__read_inp(fname)
