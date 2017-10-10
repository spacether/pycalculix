"""This module defines the Component class, which is used to make components.

Components are collections that loads are applied to.
"""
from . import base_classes
from . import geometry

class Component(base_classes.Idobj):
    """Makes a component to store loads, boundary conditions, matl, thickness.

    A component must store a list of items where the type is identical for
    all items.

    Args:
        item_list (list): list of items to store
        ctype (str): type of component. Options:

            - 'nodes': all element nodes (corner and mid-side)
            - 'n1': element corner nodes
            - 'faces': faces
            - 'elements': elements

        cname (str): component name
        write (bool): whether to write the component to the inp file
    """

    def __init__(self, item_list, ctype, cname='', write=True):
        self.items = item_list
        self.ctype = ctype
        self.name = cname+'_'+ctype
        self.write = write
        base_classes.Idobj.__init__(self)

    def __eq__(self, other):
        """ Returns bool telling if self == other."""
        if isinstance(other, Component):
            if self.items == other.items and self.ctype == other.ctype:
                return True
            return False
        else:
            return False

    def __hash__(self):
        """Returns a hash of self instance."""
        return self.id

    def set_name(self, name):
        """Set the component name to the passed name.

        Args:
          name (str): passed name, component name is changed to this
        """
        self.name = name

    def get_name(self):
        """Returns the component name created with its id number."""
        return 'C'+str(self.id)

    def get_children(self):
        """Returns a list the child ctype items of the current component."""
        res = []
        for i in self.items:
            children = getattr(i, self.ctype)
            if isinstance(children, list):
                res += children
            else:
                res += [children]
        return res

    def ccx(self):
        """Writes a component for Calculix ccx solver.

        Returns a list strings, where each item is a line to write
        to a Calculix .inp text file.

        Note:
            Only node and element components should use this.
        """
        res = []
        items_per_line = 6
        firstline = ''
        if self.ctype == 'nodes':
            # node component, displacement
            firstline = '*NSET,NSET='+self.name
        elif self.ctype == 'elements':
            # element component, materials
            firstline = '*ELSET,ELSET='+self.name
        elif self.ctype == 'faces':
            # face component, contact surface
            firstline = '*SURFACE,NAME=%s,TYPE=ELEMENT' % self.name
        res.append(firstline)
        items = self.get_children()
        if self.ctype == 'faces':
            for face in items:
                enum = face.element.id
                fnum = face.id
                res.append("%i,S%i" % (enum, fnum))
        elif self.ctype in ['nodes', 'elements']:
            grouped_items = base_classes.chunk_list(items, items_per_line)
            for group in grouped_items:
                item_ids = [str(x.id) for x in group]
                line = ', '.join(item_ids)
                if group != grouped_items[-1]:
                    line += ','
                res.append(line)
        return res

    def write_cgx(self):
        """Writes a component for Calculix cgx prerocessor.

        Returns a list strings, where each item is a line to write
        to a Calculix .fbd text file.

        Note:
            Only line and point node components should use this.
        """
        res = []

        # set the type of component
        parent = 'l'
        if (isinstance(self.items[0], geometry.Line) or
                isinstance(self.items[0], geometry.Arc)):
            parent = 'l'
        elif isinstance(self.items[0], geometry.Point):
            parent = 'p'

        # pull the list of entities
        alist = ' '.join([a.get_name() for a in self.items])

        if self.ctype == 'n1':
            # write component of line element first order nodes
            res.append('seta '+self.name+' '+parent+' '+alist)
        elif self.ctype == 'nodes':
            res.append('seta '+self.name+' '+parent+' '+alist)
            res.append('comp '+self.name+' do')
            # write component of all line element nodes
        elif self.ctype == 'f':
            # write component of line element faces
            res.append('seta '+self.name+' '+parent+' '+alist)
            res.append('comp '+self.name+' do')
            res.append('comp '+self.name+' do')
        return res

    def write_gmsh(self):
        """Writes a component for gmsh mesher.

        Returns a list strings, where each item is a line to write
        to a Calculix .geo text file.

        Note:
            Only line and point node components should use this.
        """
        res = []

        # set the type of component
        parent = 'l'
        if (isinstance(self.items[0], geometry.Line) or
                isinstance(self.items[0], geometry.Arc)):
            parent = 'l'
        elif isinstance(self.items[0], geometry.Point):
            parent = 'p'

        # pull the list of entities
        alist = ','.join([str(a.id) for a in self.items])

        line = ''
        if parent == 'l':
            line = "Physical Line('"+self.name+"') = {" + alist + '};'
        elif parent == 'p':
            line = "Physical Point('"+self.name+"') = {" + alist + '};'
        res.append(line)
        return res
