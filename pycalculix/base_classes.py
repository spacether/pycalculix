"""This module sotres base classes and functions used by others."""
from math import ceil

RESFIELDS = {}
RESFIELDS['displ'] = 'ux,uy,uz,utot'.split(',')
RESFIELDS['stress'] = 'Sx,Sy,Sz,Sxy,Syz,Szx,Seqv,S1,S2,S3'.split(',')
RESFIELDS['strain'] = 'ex,ey,ez,exy,eyz,ezx,eeqv,e1,e2,e3'.split(',')
RESFIELDS['force'] = 'fx,fy,fz'.split(',')
"""RESFIELDS is a dict which stores results fields under results type strings.

Args:
  key (str):
    'displ': displacement
    'stress': stress
    'strain': strain
    'force': force

Returns:
  list:
    'displ' = 'ux,uy,uz,utot'.split(',')
    'stress' = 'Sx,Sy,Sz,Sxy,Syz,Szx,Seqv,S1,S2,S3'.split(',')
    'strain' = 'ex,ey,ez,exy,eyz,ezx,eeqv,e1,e2,e3'.split(',')
    'force' = 'fx,fy,fz'.split(',')
"""

#FIELDTYPE is a dict that inverts the RESFIELDS dictionary mapping.
FIELDTYPE = {}
for (k, v) in RESFIELDS.items():
    for vi in v:
        FIELDTYPE[vi] = k

class Idobj(object):
    """Makes an object that stores an id number for an object.

    This is a base class for nodes, lines, areas etc.

    Atributes:
      id (int): the unique id number for the item
    """

    def __init__(self):
        self.id = -1

    def set_id(self, id):
        """Set the id number """
        self.id = id

    def __hash__(self):
        """Set the hash of the item to the id number.

        This allows one to make sets of these items.
        """
        return self.id


class Itemlist(list):
    """Makes a custom list used to store lists of non-mesh items.

    non-mesh items = Point, Line, Arc, Area, Part, Components, Loads
    This class allows us to automatically assign id numbers to items in the
    list. Thes id numbers are needed when we send geometry out for meshing.
    All instances of this class start as empty lists.
    """

    def __init__(self):
        super().__init__() # these lists start empty

    def get_ids(self):
        """Returns a list of ids. Loops through all items in self."""
        return [a.id for a in self]

    def get_next_id(self):
        """Returns the next unused id."""
        ids = self.get_ids()
        minid = 0
        if len(ids) == 0:
            return minid    # list is empty so return the minid number
        else:
            ids = sorted(ids)
            maxid = ids[-1]
            unused = list(set(list(range(minid, maxid+2))) - set(ids))
            return unused[0]    # return first available id number

    def append(self, item):
        """Adds an item to the list and sets the item's id.

        Args:
          item (Point or Line or Arc or Area or Part): item to add to the list
        """
        idnum = self.get_next_id()
        item.set_id(idnum)
        super().append(item)
        return item


class Meshlist(list):
    """Makes a custom list used to store lists of nodes and elements.

    All instances of this class start as empty lists.
    """

    def __init__(self):
        list.__init__([])

    def get_minid(self):
        """Returns the min id number in the list."""
        ids = [x.id for x in self]
        return min(ids)

    def get_maxid(self):
        """Returns the max id number in the list."""
        ids = [x.id for x in self]
        return max(ids)

    def idget(self, idnum):
        """Returns an item with the passed id number.

        Args:
          idnum (int): passed id number, we want the item that has this

        Returns:
          item or None: returns item if found, None if not found
        """
        for item in self:
            if item.id == idnum:
                return item
        #print ('ID %i was not found!' % (idnum))
        return None

    def set_minid(self, val):
        """Sets the min id to the passed val.

        Args:
          val (int): New min id number in the list
        """
        print('Re-indexing elements to min element number of: %i' % val)
        print('Old min:%i  max:%i' % (self.get_minid(), self.get_maxid()))
        minid = self.get_minid()
        offset = minid - val
        for item in self:
            item.id -= offset
        print('New min:%i  max:%i' % (self.get_minid(), self.get_maxid()))


def chunk_list(inlist, size):
    """Returns a list of lists where each list <= size length.

    Splits inlist into list of lists where each child list len <= size.

    Args:
      inlist (list): list that we want cut into smaller lists
      size (int): max length of small lists returned

    Returns:
      res (list of lists): list of list where each child list <= size length
    """
    res = []
    numlists = ceil(len(inlist)/size)
    for ind in range(numlists):
        res.append(inlist[ind*size:(ind+1)*size])
    return res
