from . import geometry

class Element(object):
    """Makes a mesh element.
    
    Supports 4 noded quads, 8 noded quads, 3 noded tris, 6 noded tris.
    
    Args:
      enum (int): element id number
      ccxtype (str): ccx element type string, 4 characters long

        Description: element_shape + elemenet_order + element_type
        element_shape = 'tri' or 'quad'
        element_order = '1' or '2'
          '1' = corner nodes only (3 or 4 noded element)
          '2' = corner nodes and midside nodes (6 or 8 noded element)
        element_shape = 'quad' or 'tri'
          'quad' = quadrangle
          'tri' = triangle
        # axisym
        'CAX6' : 'tri2axisym'
        'CAX3' : 'tri1axisym'
        'CAX8' : 'quad2axisym'
        'CAX4' : 'quad1axisym'
        # plane stress
        'CPS6' : 'tri2plstress'
        'CPS3' : 'tri1plstress'
        'CPS8' : 'quad2plstress'
        'CPS4' : 'quad1plstress'
        # plane strain
        'CPE6' : 'tri2plstrain'
        'CPE3' : 'tri1plstrain'
        'CPE8' : 'quad2plstrain'
        'CPE4' : 'quad1plstrain'
        
      nlist (list): list of Node, the nodes that define the element
    
    Attributes:
      id (int): element id numbder
      ccxtype (str): ccx element type string, 4 characters long
        ccx element type : Text description element_shape + elemenet_order + element_type
          element_shape = 'tri' or 'quad'
          element_order = '1' or '2'
            '1' = corner nodes only (3 or 4 noded element)
            '2' = corner nodes and midside nodes (6 or 8 noded element)
          element_type = specifies plane stress, plane strain, or axisymmetric
            'plstress' = plane stress
            'plstrain' = plane strain
            'axisym' = axisymmetric
        # axisym
        'CAX6' : 'tri2axisym'
        'CAX3' : 'tri1axisym'
        'CAX8' : 'quad2axisym'
        'CAX4' : 'quad1axisym'
        # plane stress
        'CPS6' : 'tri2plstress'
        'CPS3' : 'tri1plstress'
        'CPS8' : 'quad2plstress'
        'CPS4' : 'quad1plstress'
        # plane strain
        'CPE6' : 'tri2plstrain'
        'CPE3' : 'tri1plstrain'
        'CPE8' : 'quad2plstrain'
        'CPE4' : 'quad1plstrain'
      node (dict): dictionary of nodes, keys are int >= 1
      face (dict): dictionary of faces, keys are int >= 1
      center (Point): the center point of the element
    """
    def __init__(self, enum, ccxtype, nlist):
        self.id = enum
        self.ccxtype = ccxtype
        self.node = {}
        self.face = {}

        # calculate the number of faces
        fnum = len(nlist)
        if fnum > 4:
            fnum = int(fnum/2)

        # store nodes
        for (ind, node) in enumerate(nlist):
            if ind >= fnum:
                node.set_order(2)
            self.node[ind+1] = node
        
        # store faces
        n = nlist[0:fnum]+[nlist[0]]
        for i in range(fnum):
            n1 = n[i]
            n2 = n[i+1]
            f = Face(i+1, n1, n2, self)
            if len(nlist) > fnum:
                # we have second order nodes
                f.set_nmid( self.node[i+1+fnum] )
            self.face[i+1] = f
        
        # set the element center
        self.center = self.calc_center()
    
    def calc_center(self):
        """Returns the element center as a Point."""
        pts = self.get_corner_nodes()
        axials = [p[0] for p in pts]
        radials = [p[1] for p in pts]
        ax = sum(axials)/len(axials)
        rad = sum(radials)/len(radials)
        return geometry.Point(rad, ax)
            
    def get_tris(self):
        """Returns a list of triangles for plotting. Triangles are node ids.

        Triangles are closed in a counter clockwise (CCW) direction.
        """
        res = []
        numnodes = len(self.nodes()) 
        if numnodes == 3:
            # triangle, first order = 1 triangle
            res.append( [self.node[1].id, self.node[2].id, self.node[3].id] )
        elif numnodes == 4:
            # quad, first order = 2 trianges
            res.append( [self.node[1].id, self.node[2].id, self.node[3].id] )
            res.append( [self.node[1].id, self.node[3].id, self.node[4].id] )
        elif numnodes == 6:
            # tri, second order = 4 triangles
            res.append( [self.node[1].id, self.node[4].id, self.node[6].id] )
            res.append( [self.node[4].id, self.node[2].id, self.node[5].id] )
            res.append( [self.node[6].id, self.node[5].id, self.node[3].id] )
            res.append( [self.node[6].id, self.node[4].id, self.node[5].id] )
        elif numnodes == 8:
            # quad, second order = 6 triangles
            res.append( [self.node[1].id, self.node[5].id, self.node[8].id] )
            res.append( [self.node[5].id, self.node[2].id, self.node[6].id] )
            res.append( [self.node[6].id, self.node[3].id, self.node[7].id] )
            res.append( [self.node[8].id, self.node[7].id, self.node[4].id] )
            res.append( [self.node[8].id, self.node[5].id, self.node[6].id] )
            res.append( [self.node[8].id, self.node[6].id, self.node[7].id] )
        return res
    
    def get_area(self):
        """Returns the element area."""
        asum = 0.0
        for face in self.faces():
            n1 = face.nodes[0]
            n2 = face.nodes[1]
            asum += n2.y*n1.x - n2.x*n1.y
        asum = asum*0.5
        return asum
    
    def faces(self):
        """Returns a list of element faces."""
        return self.face.values()
    
    def nodes(self):
        """Returns a list of element nodes."""
        return self.node.values()
    
    def get_corner_nodes(self):
        """Returns a list of element corner nodes. [[axial, radial],...]"""
        # this is used for plotting so return [axial, radial]
        res = []
        nnum = len(self.faces())
        for ind in range(1,nnum+1):
            node = self.node[ind]
            res.append( [node.y, node.x] )
        return res
    
    def set_ccxtype(self, etype):
        """Sets the etype to the ccx element type string.

        Converts the human readable element type string into the ccx
        element type string and assigns it to this.etype
        
        Input human readable etype Examples:
        quad2axisym
        tri1plstrain

        Args:
          etype (string): input string describing needed type, human readable
            'plstress' = plane stress
            'plstrain' = plane strain
            'axisym' = axisymmetric
        """
        fnum = len(self.faces())
        shape = 'tri'
        if fnum == 4:
            shape = 'quad'
        order = '2'
        if len(self.nodes()) == fnum:
            order = '1'
        estr = shape+order+etype

        # Calculix CCX elements
        # axisymmetric
        _ccx_elements = {}
        _ccx_elements['tri2axisym'] = 'CAX6'
        _ccx_elements['tri1axisym'] = 'CAX3'
        _ccx_elements['quad2axisym'] = 'CAX8'
        _ccx_elements['quad1axisym'] = 'CAX4'
        # plane stress
        _ccx_elements['tri2plstress'] = 'CPS6'
        _ccx_elements['tri1plstress'] = 'CPS3'
        _ccx_elements['quad2plstress'] = 'CPS8'
        _ccx_elements['quad1plstress'] = 'CPS4'
        # plane strain
        _ccx_elements['tri2plstrain'] = 'CPE6'
        _ccx_elements['tri1plstrain'] = 'CPE3'
        _ccx_elements['quad2plstrain'] = 'CPE8'
        _ccx_elements['quad1plstrain'] = 'CPE4'
                
        self.ccxtype = _ccx_elements[estr]
        
    def ccx(self):
        """Returns text string defining the element for a ccx input file."""
        nids = [str(n.id) for n in self.node.values()]
        val = str(self.id)+', '+', '.join(nids)
        return val

    def get_name(self):
        """Returns the element name based on id number."""
        name = 'E%i' % (self.id)
        return name
    
    def __str__(self):
        """Returns the element name."""
        return ('E%i' % (self.id))


class Face(object):
    """Makes an element face.
    
    Args:
      fnum (int): element face number
      n1 (Node): first node
      n2 (Node): second node
      element (Element): parent element
    
    Attributes:
      id (int): element face number
      nodes (list): list of nodes [n1, n2]
      nmid (Node or None): midside node, defaults to None
      element (Element): parent element
      ext (bool): stores if face is external   
    """
    def __init__(self, fnum, n1, n2, element):
        self.id = fnum
        self.nodes = [n1,n2]
        self.nmid = None
        self.element = element
        self.ext = False

    def __hash__(self):
        """Returns a hash = self.element.id*4, Lets us make sets of faces."""
        # we assume the largest number of faces is 4, face id is 1-4, need 0-3
        hval = self.element.id*4+(self.id-1)
        return hval

    def length(self):
        """Return the length of the face."""
        p0 = geometry.Point(self.nodes[0].x,  self.nodes[0].y)
        p1 = geometry.Point(self.nodes[1].x,  self.nodes[1].y)
        vect = p1 - p0
        return vect.length()

    def get_mnorm(self):
        """Returns a list: [midpoint, noarmal_vector]

        Returns:
          midpoint (Point): face midpoint
          norm_vect (Point): face normal vector, it is a unit vector
        """
        # return midpt, normal
        p0 = geometry.Point(self.nodes[0].x,  self.nodes[0].y)
        p1 = geometry.Point(self.nodes[1].x,  self.nodes[1].y)
        midpoint = p0 + p1
        midpoint = midpoint*0.5
        norm_vect = p1 - p0
        norm_vect.rot_ccw_deg(90)
        norm_vect.make_unit()
        return [midpoint, norm_vect]

    def set_ext(self):
        """Sets the ext boolean flag to True. Face is external."""
        self.ext = True

    def set_nmid(self, nmid):
        """Sets the face midside node, nmid.
        
        Args:
          nmid (Point): passed node to set as the face midside node
        """
        self.nmid = nmid

    def __eq__(self, other):
        """Compare this face to other face. True if node sets equal.
        
        Args:
          other (Face): the other face we are comparing
        
        Returns:
          bool: True if face node sets are equal, False otherwise.
        """
        my_nodes = set([self.nodes[0].id,self.nodes[1].id])
        other_nodes = set([other.nodes[0].id,other.nodes[1].id])
        if my_nodes == other_nodes:
            return True
        else:
            return False

    def __str__(self):
        """Returns string listing object type, id number and nodes."""
        n1 = self.nodes[0].id
        n2 = self.nodes[1].id
        fid = self.id
        eid = self.element.id
        s = 'Face %i on element %i with nodes [%i,%i]' % (fid, eid, n1, n2)
        return s


class Node(object):
    """Makes a mesh node.
    
    Args:
      node_num (int): node number
      x (float): x-coordinate
      y (float): y-coordinate
      z (float): z-coordinate
    
    Attributes:
      id (int): node number
      x (float): x-coordinate
      y (float): y-coordinate
      z (float): z-coordinate
      order (int): 1 or 2, 1 if it is a corner node, 2 if it is a midside node
      elements (set): a set of elements that contain this node    
    """
    def __init__(self, node_num, x, y, z):
        self.id = node_num
        self.x = x
        self.y = y
        self.z = z
        self.order = 1
        self.elements = set()

    def set_order(self, order):
        """Sets the node order, must be 1 or 2.
        
        Args:
          order (int): node order, 1 is corner node, 2 is midside node
        """
        self.order = order
    
    def add_element(self, e):
        """Adds an element to this node's set of elements.
        
        Args:
          e (Element): element to add to set self.elements
        """
        self.elements.add(e)

    def __hash__(self):
        """Returns the node id as the hash."""
        return self.id

    def __eq__(self, other):
        """Returns boolean of id equality to other Node.
        
        Args:
          other (Node): the node to compare to this node.
        
        Returns:
          boolean: True if ids are equal, False if they are not
        """
        if isinstance(other, Node):
            return self.id == other.id
        else:
            return False
        
    def ccx(self):
        """Returns a string defining the node for a ccx inp file."""
        val = '%i, %f, %f, %f' % (self.id, self.x, self.y, self.z)
        return val
    
    def get_name(self):
        """Returns the string node name."""
        name = 'N%i' % (self.id)
        return name

    def __str__(self):
        """Returns string listing object type, id number and (x,y) coords."""
        val = 'Node, id %i, (x,y)=(%f,%f)' % (self.id, self.x,self.y)
        return val
    