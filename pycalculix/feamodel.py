import matplotlib.pyplot as plt
from matplotlib.patches import Polygon  # needed for plotting elements
from matplotlib.collections import PatchCollection  # element plotting
import matplotlib.colors as colors
import matplotlib.cm as cmx
import subprocess # used to launch meshers cgx and gmsh
import os # used to delete files written by cgx
from numpy import linspace

from . import environment
from . import base_classes
from . import geometry
from . import material
from . import components
from . import loads
from . import mesh
from . import part
from . import model

# element colors, 0-1, 0=black, 1=whate
_ecolor = '.4'
_fcolor = '.9'

class FeaModel(object):
    def __init__(self, fname, ccx=None, cgx=None, gmsh=None):
        self.fname = fname
        self.points = base_classes.Itemlist()
        self.lines = base_classes.Itemlist()
        self.areas = base_classes.Itemlist()
        self.parts = base_classes.Itemlist()
        self.matls = base_classes.Itemlist()
        self.components = base_classes.Itemlist()
        self.loads = {}    # store loads by time
        self.models = base_classes.Itemlist()
        self.nodes = base_classes.Meshlist()
        self.elements = base_classes.Meshlist()
        self.faces = []
        self.sel = {}   # this selected set is what we'll plot
        self.time = 1.0 # 0 is model set-up, 1 is first step, etc
        self.units = {}
        self.select()
        
        # fix the paths to the needed programs, ccx, cgx, and gmsh
        if ccx != None:
            environment.CCX = ccx
        if cgx != None:
            environment.CGX = cgx
        if gmsh != None:
            environment.GMSH = gmsh
        
    def select(self, items = 'all', byfaces=True):
        # set the selected item to the passed item or list of elements
        # if byfaces is set to true, only elements with faces selected will be
        # shown

        allsel = False
        if isinstance(items, str):
            if items == 'all':
                allsel = True
                self.sel['nodes'] = self.nodes
                self.sel['elements'] = self.elements
                self.sel['points'] = self.points
                self.sel['lines'] = self.lines
                self.sel['areas'] = self.areas
                print('Selected all')                

        items = self.listify(items)
        self.sel['selected'] = items

        if allsel == False:
            # select children
            mystr = [str(a) for a in items]
            mystr = '+'.join(mystr)
            print('Selected: '+mystr)
            self.sel_children(byfaces)

    def sel_children(self, byfaces=True):
        # this selects the children of the selected items
        types = ['nodes','elements','faces','points','lines','areas']
        # zero out the selection dict
        for t in types:
            self.sel[t] = set()
        for item in self.sel['selected']:
            if isinstance(item, mesh.Node):
                self.sel['nodes'].add(item)
                
            elif isinstance(item, mesh.Element):
                self.sel['elements'].add(item)
                child_nodes = item.nodes()
                self.sel['nodes'].update( child_nodes )
                
            elif isinstance(item, mesh.Face):
                self.sel['faces'].add(item)
                
            elif isinstance(item, geometry.Point):
                self.sel['points'].add(item)
                for n in item.nodes:
                    self.sel['elements'].update(n.elements)
                for e in self.sel['elements']:
                    self.sel['nodes'].update( e.nodes() )
                    
            elif (isinstance(item, geometry.Line) or 
                  isinstance(item, geometry.Arc)):
                self.sel['points'].update(item._points_)
                self.sel['lines'].add(item)
                self.sel['faces'].update(item.faces)
                if byfaces:
                    # add elements from face only
                    for f in item.faces:
                        self.sel['elements'].add(f.element)
                        self.sel['nodes'].update( f.element.nodes() )
                else:
                    # add elements connected to selected nodes, add their nodes 
                    self.sel['nodes'].update(item.nodes)
                    for n in item.nodes:
                        self.sel['elements'].update(n.elements)
                    for e in self.sel['elements']:
                        self.sel['nodes'].update( e.nodes() )
                        
            elif isinstance(item, geometry.Area):
                pts = item.get_points()
                self.sel['points'].update(pts)
                self.sel['lines'].update(item.lines)
                for l in item.lines:
                    self.sel['faces'].update(l.faces)
                self.sel['nodes'].update(item.nodes)
                self.sel['elements'].update(item.elements)
                
            elif isinstance(item, part.Part):
                pts = item.get_points()
                lines = item.get_lines()
                self.sel['points'].update(pts)
                self.sel['lines'].update(lines)
                for l in lines:
                    self.sel['faces'].update(l.faces)
                self.sel['nodes'].update(item.nodes)
                self.sel['elements'].update(item.elements)

    def set_units(self, dist='m', temp='K'):
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
        self.units = d
        
    def get_units(self, *args):
        # returs units for the passed argumnets
        
        res = []
        for a in args:
            mystr = ''
            if a in self.units:
                # check if the requested item is in the units dict
                mystr = ' ('+ self.units[a] + ')'
            elif a in base_classes.FIELDTYPE:
                # return units on results dict item, ex: stress, strain etc
                ftype = base_classes.FIELDTYPE[a]
                if ftype in self.units:
                    mystr = ' ('+ self.units[ftype] + ')'
            res.append(mystr)
        return res
        
    def set_time(self, time):
        # sets the time in the feamodel (preprocessor)
        self.time = time
        
    def get_item(self, item):
        # return item from part given string identifying item
        if item[0] == 'P':
            # get point
            items = self.points
            num = int(item[1:])
            res = [a for a in items if a.id == num]
            return res[0]
        elif item[0] == 'L':
            # get line
            items = self.lines
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

    def MatlMaker(self, name):
        mat = material.Material(name)
        self.matls.append(mat)
        return mat

    def PartMaker(self):
        p = part.Part(self)
        self.parts.append(p)
        return p

    def ModelMaker(self, parts, mtype):
        m = model.Model(self, parts, mtype)
        self.models.append(m)
        return m    

    def plot_elements(self, fname='', items=None, display=True, enum=False, nnum=False):
        # plots the elements in the part

        nodes = []
        elements = []
        if items == None:
            # use selected set
            items = self.sel['selected']
            nodes = self.sel['nodes']
            elements = self.sel['elements']
        else:
            # used passed items
            items = self.listify(items)
            for i in items:
                nodes += i.nodes
                elements += i.elements
        
        if len(elements) > 0:
                
            radials = [n.x for n in nodes]
            axials = [n.y for n in nodes]
            vadder = (max(radials) - min(radials))/5
            hadder = (max(axials) - min(axials))/5
    
            (vmax,vmin) = (max(radials)+vadder, min(radials)-vadder)
            (hmax,hmin) = (max(axials)+hadder, min(axials)-hadder)
        
            # plotting elements
            fig, ax = plt.subplots()            
            patches = []
            for e in elements:
                xycoords = e.get_corner_nodes()
                polygon = Polygon(xycoords, closed=True)
                patches.append(polygon)
            p = PatchCollection(patches, facecolors=_fcolor, edgecolors=_ecolor)
            ax.add_collection(p)
            
            # add element numbers
            if enum:
                for e in elements:
                    [ename, pax, prad] = [e.get_name(), e.center.y, e.center.x]
                    ax.annotate(ename, (pax,prad))

            # add element numbers
            if nnum:
                for n in nodes:
                    [nname, pax, prad] = [n.get_name(), n.y, n.x]
                    ax.annotate(nname, (pax,prad))

            
            # set units
            [d_unit] = self.get_units('dist')            
            
            iname = self.get_cname(items)
            plt.title(iname+' elements')
            plt.xlabel('axial, y'+d_unit)
            plt.ylabel('radial, x'+d_unit)
            plt.axis('scaled')
            plt.xlim(hmin, hmax)
            plt.ylim(vmin, vmax)
            
            if fname != '':
                # save the image
                fname += '.png'
                if environment.DPI != None:
                    plt.savefig(fname, dpi=environment.DPI, bbox_inches='tight')
                else:
                    plt.savefig(fname, bbox_inches='tight')
            
            if display:
                plt.tight_layout()
                plt.show()        

            # remove all figures
            plt.close()

        else:
            # part has not been meshed yet
            res = 'No elements exist! '
            res += 'Try meshing your parts. ' 
            print(res)

    def plot_pressures(self, fname='', items = None, display=True):
        # plot the pressures on elements

        nodes = []
        elements = []
        if items == None:
            # use selected set
            items = self.sel['selected']
            nodes = self.sel['nodes']
            elements = self.sel['elements']
        else:
            # use passed item
            items = self.listify(items)
            for i in items:
                nodes += i.nodes
                elements += i.elements
        
        # plot all elements and store length, length determins min pressure arrow
        if len(elements) > 0:
            
            # get axials and radials for bounds
            radials = [p.x for p in nodes]
            axials = [p.y for p in nodes]
            
            # plotting elements
            fig = plt.figure()
            ax = fig.add_subplot(111,aspect='equal')

            patches = []
            face_len = []
            for e in elements:
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
            for load in self.loads[self.time]:
                if load.ltype in  ['press','press_fluid']:
                    plist += load.get_list()
            pressures = [pval for [face,pval] in plist]

            # check min and max bounds
            pmin = min(pressures)
            pmax = max(pressures)
            mult = 1
            if pmin == 0:
                mult = 1
            else:
                mult = face_len/abs(pmin)    # mult to go from pressure to length
            
            # make tick list for later plot, and color map
            tick_list = [pmin]  # default to plot one val
            cmap = colors.ListedColormap(['b', 'b']) # default to plot one val
            if pmax != pmin:
                # we have a range of values we're plotting
                tick_list = linspace(pmin,pmax,8)
                cmap = plt.cm.jet
                        
            # set color contours for arrows
            cNorm  = colors.Normalize(vmin=pmin, vmax=pmax)
            scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=cmap)            
            scalarMap._A = [] # need to set this for it to work
            
            # make arrows
            for [face, pval] in plist:
                [p1, unit] = face.get_mnorm()
                pdelta = None
                if pmin == 0:
                    pdelta = unit*face_len
                else:
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
            [d_unit, p_unit, t_unit] = self.get_units('dist', 'stress', 'time')

            # set plot axes
            iname = self.get_cname(items)
            tstr = '%s pressures %s\nTime=%f%s' % (iname, p_unit, self.time,
                                                 t_unit)
            plt.title(tstr)
            plt.xlabel('axial, y'+d_unit)
            plt.ylabel('radial, x'+d_unit)

            # set min and max vertical and axial limits
            plt.xlim(hmin, hmax)
            plt.ylim(vmin, vmax)            
            
            # set the colorbar
            plt.colorbar(scalarMap, orientation='vertical', ticks=tick_list)
            
            if fname != '':
                # save the image
                fname += '.png'
                if environment.DPI != None:
                    plt.savefig(fname, dpi=environment.DPI, bbox_inches='tight')
                else:
                    plt.savefig(fname, bbox_inches='tight')
            
            if display:
                plt.tight_layout()
                plt.show()        

            # remove all figures
            plt.close()

        else:
            # part has not been meshed yet
            res = 'Part: %s does not have any elements! ' % (self)
            res += 'Try meshing it with model.mesh(1)' 
            print(res)

    def plot_constraints(self, fname='', items = None, display=True):
        # this plots the constraints on the passed part or selected items

        nodes = []
        elements = []
        if items == None:
            # use selected set
            items = self.sel['selected']
            nodes = self.sel['nodes']
            elements = self.sel['elements']
        else:
            # use passed item
            items = self.listify(items)
            for i in items:
                nodes += i.nodes
                elements += i.elements
        
        # plot all elements
        if len(elements) > 0:
            
            # get axials and radials for bounds
            radials = [p.x for p in nodes]
            axials = [p.y for p in nodes]
            
            # plotting elements
            fig = plt.figure()
            ax = fig.add_subplot(111,aspect='equal')

            patches = []
            face_len = []            
            for e in elements:
                face_len.append( e.face[1].length() )                
                xycoords = e.get_corner_nodes()
                polygon = Polygon(xycoords, closed=True)
                patches.append(polygon)
            p = PatchCollection(patches, edgecolors=_ecolor, facecolors=_fcolor)
            ax.add_collection(p)

            # average face length is min arrow length
            face_len = sum(face_len)/len(face_len)

            # store displacements we'll plot
            # this is a list of lists [node, dict ux,uy,uz]
            ulist = []
            vals = []
            for load in self.loads[self.time]:
                if load.ltype in  ['ux','uy','uz']:
                    ulist += load.get_list()
                    # still need to filter out items not in our selection set
                    for udisp in ulist:
                        vals += list( udisp[1].values() )

            # check min and max bounds
            pmin = min(vals)
            pmax = max(vals)
            
            # make tick list for later plot, and color map
            tick_list = [pmin]  # default to plot one val
            cmap = colors.ListedColormap(['b', 'b']) # default to plot one val
            if pmax != pmin:
                # we have a range of values we're plotting
                tick_list = linspace(pmin,pmax,8)
                cmap = plt.cm.jet
                        
            # set color contours for arrows
            cNorm  = colors.Normalize(vmin=pmin, vmax=pmax)
            scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=cmap)            
            scalarMap._A = [] # need to set this for it to work
            
            # make arrows for displacements
            pvect= {'ux':geometry.Point(1,0,0),
                    'uy':geometry.Point(0,1,0),
                    'uz':geometry.Point(0,0,1)}
            for [node, udict] in ulist:
                for (k, v) in udict.items():
                    delta = pvect[k]
                    pstart = geometry.Point(node.x, node.y, node.z)
                    pend = pstart + delta*v
                    hw = face_len*0.4
                    hl = hw
                    if v == 0:
                        pend = geometry.Point(node.x, node.y, node.z)
                        pstart = pend - delta*(hl*2)
                        #radials.append(pstart.x)
                        #axials.append(pstart.y)
                    else:
                        #radials.append(pend.x)
                        #axials.append(pend.y)
                        pass
                
                    delta = pend - pstart
                    colorVal = scalarMap.to_rgba(v)
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
            [d_unit, t_unit] = self.get_units('dist', 'time')

            # set plot axes
            iname = self.get_cname(items)
            tstr = '%s constraints %s\nTime=%f%s' % (iname, d_unit, self.time,
                                                 t_unit)
            plt.title(tstr)
            plt.xlabel('axial, y'+d_unit)
            plt.ylabel('radial, x'+d_unit)

            # set min and max vertical and axial limits
            plt.xlim(hmin, hmax)
            plt.ylim(vmin, vmax)            
            
            # set the colorbar
            plt.colorbar(scalarMap, orientation='vertical', ticks=tick_list)
            
            if fname != '':
                # save the image
                fname += '.png'
                if environment.DPI != None:
                    plt.savefig(fname, dpi=environment.DPI, bbox_inches='tight')
                else:
                    plt.savefig(fname, bbox_inches='tight')
            
            if display:
                plt.tight_layout()
                plt.show()        

            # remove all figures
            plt.close()

        else:
            # part has not been meshed yet
            res = 'Part: %s does not have any elements! ' % (self)
            res += 'Try meshing it with model.mesh(1)' 
            print(res)

    def plot_geometry(self, fname='', items = None, display=True):
        # this method plots the part
        # check out: http://nickcharlton.net/posts/drawing-animating-shapes-matplotlib.html

        points = []
        lines = []
        areas = []
        if items == None:
            # use selected set
            items = self.sel['selected']
            points = self.sel['points']
            lines = self.sel['lines']
            areas = self.sel['areas']
        else:
            items = self.listify(items)
            if isinstance(items[0], part.Part):
                for i in items:
                    points += i.get_points()
                    lines += i.get_lines()
                    areas += i.areas
            elif isinstance(items[0], geometry.Area):
                for i in items:
                    points += i.get_points()
                    lines += i.lines
                    areas += i
            elif isinstance(items[0], geometry.Line):
                for i in items:
                    points += i._points_
                    lines += i                
            elif isinstance(items[0], geometry.Point):
                for i in items:
                    points += i
        
        nax=[pt.y for pt in points]
        nrad=[pt.x for pt in points]
        n_id=[pt.get_name() for pt in points]
        
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
        for l in lines:
            l.plot(ax)
        
        #plot area ids
        for a in self.areas:
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
        [d_unit] = self.get_units('dist')            
        
        # show plot
        iname = self.get_cname(items)
        plt.title(iname+' geometry')
        plt.xlabel('axial, y'+d_unit)
        plt.ylabel('radial, x'+d_unit)
        ax.set_aspect('equal')

        if fname != '':
            # save the image
            fname += '.png'
            if environment.DPI != None:
                plt.savefig(fname, dpi=environment.DPI, bbox_inches='tight')
            else:
                plt.savefig(fname, bbox_inches='tight')
        
        if display:
            plt.tight_layout()
            plt.show()
        
        # remove all figures
        plt.close()
    
    def listify(self, items):
        # converts item into a list if it's not one
        if not isinstance(items, list):
            items = [items]
        return items
    
    def get_cname(self, items):
        # this returns a component name prefix, for labeling lists of items
        cname = ''
        if len(items) == 1:
            if items[0] == 'all':
                cname = 'all'
            else:
                cname = items[0].get_name()
        else:
            cname = items[0].get_name()+'-'+items[-1].get_name()
        return cname

    def add_load(self, load, time):
        # add load to feamodel and part
        if time in self.loads:
            self.loads[time].append(load)
        else:
            self.loads[time] = [load]

    def set_gravity(self, grav, items):
        # applies gravity to the items
        items = self.listify(items)
        ctype = 'elements'        
        cname = self.get_cname(items)
        comp = components.Component(items, ctype, cname)
        
        if comp not in self.components:
            self.components.append(comp)
        
        ltype = 'gravity'
        time = self.time
        load = loads.Load(ltype, comp, grav)
        
        # add load to feamodel
        self.add_load(load, time)

    def set_rpm(self, rpm, items):
        # applies rpm to the items
        items = self.listify(items)
        ctype = 'elements'
        cname = self.get_cname(items)
        comp = components.Component(items, ctype, cname)
        
        if comp not in self.components:
            self.components.append(comp)
        
        ltype = 'rpm'
        time = self.time
        load = loads.Load(ltype, comp, rpm)
        
        # add load to feamodel
        self.add_load(load, time)

    def set_radps(self, rpm, items):
        # applies radps to the items
        items = self.listify(items)
        ctype = 'elements'
        cname = self.get_cname(items)
        comp = components.Component(items, ctype, cname)
        
        if comp not in self.components:
            self.components.append(comp)
        
        ltype = 'radps'
        time = self.time
        load = loads.Load(ltype, comp, rpm)
        
        # add load to feamodel
        self.add_load(load, time)
        
    def set_fluid_press(self, str_items, rho, g, xo, po):
        # this sets pressure from water + atmosphere
        # p = po + rho*g*h, it assumes pressure increases in -x

        items = str_items
        # convert string into item(self)
        if isinstance(str_items, str):
            items = self.get_item(str_items)
        items = self.listify(items)

        ctype = 'faces'
        cname  = self.get_cname(items)
        comp = components.Component(items, ctype, cname)
        
        ltype = 'press_fluid'
        mult = rho*g
        load = loads.Load_linear(ltype, comp, po, mult, xo)
        
        # add load to feamodel and part
        self.add_load(load, self.time)            

    def set_load(self, ltype, str_items, lval, ldir=None):
        # applies a load or pressure to lines
        # ltype = 'press' or 'force'
        # str_items is a string defining a line, point or side of the part
        #    or an item or a list of items
        # ldir = 'x' or 'y' or 'z'
        # positive pressure sign is compressing the surface
        # this is used to apply them to the child nodes/faces
        
        # make component if it doesn't exist
        
        # convert string into item(self)
        items = str_items
        # convert string into item(self)
        if isinstance(str_items, str):
            items = self.get_item(str_items)
        items = self.listify(items)
        
        ctype = 'nodes'
        if ltype == 'press':
            ctype = 'faces'
        cname  = self.get_cname(items)
        comp = components.Component(items, ctype, cname)
        
        # check if component exists, and if so use it, need to write it here
        self.components.append(comp)

        if ltype == 'force':
            ltype = 'f'+ldir # for example fx

        # add load to feamodel
        load = loads.Load(ltype, comp, lval)
        self.add_load(load, self.time)

    def set_constr(self, ltype, line_list, ldir, lval=0.0):
        # this sets constraints on lines
        # ltype = 'fix' or 'displ'
        # ldir = 'x' or 'y' or 'z'
        ctype = 'nodes'

        # convert string into item(s)
        if isinstance(line_list, str):
            line_list = self.get_item(line_list)
        line_list = self.listify(line_list)

        cname  = self.get_cname(line_list)
        comp = components.Component(line_list, ctype, cname)

        # check if component exists, and if so use it, need to write it here
        self.components.append(comp)

        ltype = 'u'+ldir # for example ux
        load = loads.Load(ltype, comp, lval)

        # add load to feamodel and part
        self.add_load(load, self.time)

    def set_eshape(self, eshape='quad', eorder=2):
        # sets eleemnt properties and thickness if applicable
        self.eshape = eshape # quad or tri
        self.eorder = eorder # 1 or 2

    def set_etype(self, items, etype='plstress', thick=None):
        # sets the element type on part areas, or a list of areas

        # convert string into item(s)
        if isinstance(items, str):
            items = self.get_item(items)

        # convert items into areas
        items = self.listify(items)
        cname = self.get_cname(items)
        if isinstance(items[0], part.Part):
            tmp = []
            for ind, p in enumerate(items):
                tmp += p.areas
            items = tmp

        # set the element types on the areas, this is used to fix
        # elements when importing them from the inp file
        for area in items:
            area.set_etype(etype)

        # set a thickness component if needed
        if etype != 'axisym' and thick != None:
            
            # make component for nodal thickness
            ctype = 'nodes'
            comp = next((c for c in self.components if c.name == cname+'_nodes'), None)
            if comp == None:
                comp = components.Component(items, ctype, cname)
                self.components.append(comp)            

            # check existing thickness components and remove new area from them
            for time in self.loads:
                for load in self.loads[time]:
                    if load.ltype == 'thickness':
                        c = load.comp
                        areas = c.items
                        for newarea in items:
                            if newarea in areas:
                                areas.remove(newarea)
                        # relabel modified component
                        c.name = self.get_cname(areas)+'_nodes'                        
                        
            ltype = 'thickness'
            time = 0.0
            load = loads.Load(ltype, comp, thick)
            self.add_load(load, time)

    def set_matl(self, matl, items):
        # sets the matl on an item, part, or area + makes a component to apply it
        items = self.listify(items)
        cname = self.get_cname(items)
        if isinstance(items[0], part.Part):
            tmp = []
            for ind, p in enumerate(items):
                tmp += p.areas
            items = tmp

        # set the matls on the areas
        for area in items:
            pass
            #area.set_matl(matl)

        # make component
        ctype = 'elements'
        comp = next((c for c in self.components if c.name == cname+'_elements'), None)
        if comp == None:
            comp = components.Component(items, ctype, cname)
            self.components.append(comp)            

        # check existing matl components and remove new area from them
        for time in self.loads:
            for load in self.loads[time]:
                if load.ltype == 'matl':
                    c = load.comp
                    areas = c.items
                    for newarea in items:
                        if newarea in areas:
                            areas.remove(newarea)
                    # relabel modified component
                    c.name = self.get_cname(areas)+'_elements'                        
        
        # manually set all elments to the correct matl if they already exist
        for area in items:
            if hasattr(area, 'elements'):
                for e in area.elements:
                    pass
        
        ltype = 'matl'
        time = 0.0
        load = loads.Load(ltype, comp, matl)
        self.add_load(load, time)

    def mesh(self, fineness, mesher='cgx'):
        # this initiates meshing of all parts
        self.mesher = mesher
                
        if mesher == 'gmsh':
            self.mesh_gmsh(fineness)
        elif mesher == 'cgx':
            self.mesh_cgx(fineness)

    def mesh_gmsh(self, fine):
        # this meshes the model
        
        geo = []

        # write all points
        for pt in self.points:
            linestr = 'Point(%i) = {%f, %f, %f};' % (pt.id, pt.x, pt.y, 0.0)
            geo.append(linestr)
        
        # write all lines
        for line in self.lines:
            ln = line.id
            p1 = line.pt(0).id
            p2 = line.pt(1).id
            linestr = ''
            if isinstance(line, geometry.Arc):
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
                if self.eshape == 'quad':
                    ndiv = line.ediv/2+1
                    esize = esize*2
                    # this is needed because quad recombine
                    # splits 1 element into 2
                linestr = 'Transfinite Line{%i} = %i;' % (ln, ndiv)
                print('LINE ELEMENT SIZE: %f, MAKES %i ELEMENTS' % (line.length()/line.ediv, line.ediv))
                geo.append(linestr)
                geo.append('Characteristic Length {%i,%i} = %f;' % (p1, p2, esize))            

        # write all areas
        for area in self.areas:
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
        for p in self.parts:
            # make components for each part
            line = "Physical Surface('%s') = " % (p.get_name())
            area_ids = []
            for area in p.areas:
                if area.closed:
                    area_ids.append(str(area.id))
            line = line + '{' + ','.join(area_ids) + '};'
            geo.append(line)
                    
        # write all line componenets so we can get nodes out
        for L in self.lines:
            line = "Physical Line('%s') = {%i};" % (L.get_name(), L.id)
            geo.append(line)
        
        # write node componenets
        # node list is not produced by gmsh
        for pt in self.points:
            linestr = "Physical Point('%s') = {%i};" % (pt.get_name(), pt.id)
            geo.append(linestr)                    

        # set the meshing options
        geo.append('Mesh.CharacteristicLengthFactor = '+str(fine)+'; //mesh fineness')
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
        f = open(fname,'w')
        for line in geo:
            #print (line)
            f.write(line+'\n')
        f.close()
        print ('File: '+ fname + ' was written')
        
        # run file in bg mode, -2 is 2d mesh
        astr = "%s %s -2 -o %s" % (environment.GMSH, fname, fout)
        print(astr)
        subprocess.call(astr, shell=True)
        print ('File: '+ fout + ' was written')
        print('Meshing done!')
        
        '''
        # fix element type
        # define the element types, and tells the program to mesh
        estr = self.eshape+str(self.eorder)+self.parts[0].areas[0].etype
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
        '''

        # write gmsh msh file
        astr = "%s %s -2 -o %s" % (environment.GMSH, fname, self.fname+'.msh')
        subprocess.call(astr, shell=True)
        print ('File: '+ self.fname+ '.msh was written')

        # read in the calculix mesh
        self.read_inp(self.fname+'.inp')

    def read_inp(self, fname):
        # this function reads in the calculix input file
        # stores all nodes, elements, and faces
        # and assigns all element and node sets correctly

        f = open(fname,'r')
        mode = None
        set_name = None
        set_type = None
        
        items = [] # holder for nodes or elements in nsets or esets
        N = base_classes.Meshlist() # store nodes        
        E = base_classes.Meshlist() # store elements, allows renumbering before putting int model
        F = [] # store faces
        sets = {'E':{},'N':{}} # store sets
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
                    for n in nlist:
                        n.add_element(e)
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
        self.elements = E
        self.faces = F
        
        # remove arc center ndoes from imported node set
        for line in self.lines:
            if isinstance(line, geometry.Arc):                
                pt = line.actr
                ndist = []
                for n in N:
                    p_tmp = geometry.Point(n.x, n.y)
                    p_tmp = pt - p_tmp
                    dist = p_tmp.length()
                    ndist.append( {'dist':dist,'node':n} )
                # sort the list by dist, sorts low to high
                ndist = sorted(ndist, key=lambda k: k['dist'])
                match_node = ndist[0]['node']
                N.remove( match_node )
        
        self.nodes = N
                
        for p in self.parts:
            # assign part element and node sets
            pname = p.get_name()
            p.elements = sets['E'][pname]
            p.nodes = sets['N'][pname]
            
            # assign all nodes and elements to areas, fix the element types
            for area in self.areas:
                aname = area.get_name()
                area.elements = sets['E'][aname]
                area.nodes = sets['N'][aname]
                area.set_child_ccxtypes() #paint element types on elements
            
            # assign the child nodes to points
            pts = p.get_points()
            for pt in pts:
                ndist = []
                for n in p.nodes:
                    p_tmp = geometry.Point(n.x, n.y)
                    p_tmp = pt - p_tmp
                    dist = p_tmp.length()
                    ndist.append( {'dist':dist,'node':n} )
                # sort the list by dist, sorts low to high
                ndist = sorted(ndist, key=lambda k: k['dist'])
                pt.nodes = [ndist[0]['node']]
                #print('Point %s = node %s' % (pt, pt.nodes))
            
            # assign the nodes and n1 and faces to lines
            lines = p.get_lines()
            for line in lines:
                lname = line.get_name()
                nodes = sets['N'][lname]
                n1 = [n for n in nodes if n.order == 1]
                faces = []
                for face in external:
                    if set(face.nodes).issubset(set(nodes)):
                        faces.append(face)
                line.nodes = nodes
                line.n1 = n1
                line.faces = faces
                
        print('Done reading Calculix/Abaqus .inp file')
        self.select() # select all, adds nodes and elements to the selected set

    def mesh_cgx(self, fine):
        # meshes the parts using calculix preprocessor
        fbd = []
        comps = []
        cfiles = []
        
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
        
        num = 1.0/fine
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
            ln = line.get_name()
            p1 = line.pt(0).get_name()
            p2 = line.pt(1).get_name()
            linestr = ''
            if isinstance(line, geometry.Arc):
                # line is arc
                pc = line.actr.get_name()
                linestr = 'line %s %s %s %s' % (ln, p1, p2, pc)
            else:
                # straight line
                linestr = 'line %s %s %s' % (ln, p1, p2)
            # set division if we have it
            if hasattr(line, 'ediv'):
                for p in self.parts:
                    if line in p.get_lines():
                        break
                ndiv = p.eorder*line.ediv 
                linestr += ' '+str(int(ndiv))
            fbd.append(linestr)
            L = 'seta %s l %s' % (ln, ln)
            comps.append(L)
            cfiles.append(ln)
        
        # write all areas
        for area in self.areas:
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
        for p in self.parts:
            # make components for each part, syntax
            # seta P0 s s0 s1
            line = 'seta %s s ' % (p.get_name(),)
            cfiles.append(p.get_name())
            area_ids = []
            for area in p.areas:
                if area.closed:
                    area_ids.append(area.get_name())
            line = line + ' '.join(area_ids)
            fbd.append(line)
            #lines that mesh the part
            part_comp = p.get_name()
            
            # define the element types, and tells the program to mesh
            estr = p.eshape+str(p.eorder)+p.etype
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
        fname = self.fname+'.fbd'
        f = open(fname,'w')
        for line in fbd:
            #print (line)
            f.write(line+'\n')
        f.close()
        print ('File: '+ fname + ' was written')
        
        # run file in bg mode
        p = subprocess.call("%s -bg %s" % (environment.CGX, fname), shell=True)
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
        fname = self.fname+'.inp'
        f = open(fname,'w')
        for line in inp:
            #print (line)
            f.write(line+'\n')
        f.close()
        print ('File: '+ fname + ' was written')

        # read in the calculix mesh
        self.read_inp(fname)
