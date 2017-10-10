"""This module stores load classes.
"""
from math import pi

class ConstLoad(object):
    """Makes a load. Many types are supported.

    Possible load types:
        - forces
        - displacements
        - thickness (on nodes)
        - matl (on elements)
        - pressure (or stress on faces)
        - gravity (on elements)
        - rotation about an axis (on elements)

    Args:
        ltype (string): load type:

            - 'fx','fy','fz': force on each axis
            - 'ux','uy','uz': displacement on each axis
            - 'nodal_thickness': thickness on element nodes
            - 'matl': matl on elements
            - 'press': pressure, + is tension, - is compresion
            - 'gravity': gravity in x axis direction, - goes towards y axis
            - 'radps': radians per second rotation
            - 'rpm': rotations per minute rotation

        comp (Component): component that the load is applied to
        val (double or Matl): value of the load
            Matl is needed when setting material loads
            otherwise a double value is used to describe the load
    """

    def __init__(self, ltype, comp, val):
        self.ltype = ltype
        self.comp = comp
        self.val = val

    def get_list(self):
        """Returns a list of [item, val] for plotting."""
        res = []
        if self.ltype == 'press':
            # pressure
            faces = self.comp.get_children()
            for face in faces:
                res.append([face, self.val])
        elif self.ltype in ['ux', 'uy', 'uz']:
            nodes = self.comp.get_children()
            load_dict = {self.ltype:self.val}
            for node in nodes:
                res.append([node, load_dict])
        return res

    def ccx(self):
        """Writes a load for Calculix ccx solver.

        Returns a list strings, where each item is a line to write
        to a Calculix .inp text file.
        """
        res = []
        cname = self.comp.name
        if self.ltype == 'press':
            # pressure
            res.append('** Pressure on component: '+cname)
            res.append('*DLOAD')
            faces = self.comp.get_children()
            for face in faces:
                myline = '%i,P%i,%f' % (face.element.id, face.id, self.val)
                res.append(myline)
        elif self.ltype == 'gravity':
            # gravity
            res.append('** Gravity on component: '+cname)
            res.append('*DLOAD')
            res.append('%s,GRAV,%f,-1.,0.,0.' % (cname, self.val))
        elif self.ltype[0] == 'u':
            # displacement
            axis = self.ltype[1]
            anum = {'x':1, 'y':2, 'z':3}
            axis = anum[axis]
            res.append('*BOUNDARY')
            res.append('%s,%i,%i,%f' % (cname, axis, axis, self.val))
        elif self.ltype[0] == 'f':
            # force
            res.append('** Force on component: '+cname)
            nodes = len(self.comp.items)
            fval = self.val/nodes
            axis = self.ltype[1]
            anum = {'x':1, 'y':2, 'z':3}
            axis = anum[axis]
            res.append('*CLOAD')
            res.append('%s,%i,%f' % (cname, axis, fval))
        elif self.ltype == 'nodal_thickness':
            res.append('*NODAL THICKNESS')
            res.append('%s, %f' % (cname, self.val))
        elif self.ltype == 'radps':
            radsq = self.val**2 # square the speed in radians per sec
            res.append('** Rotation force on component: '+cname)
            res.append('** Speed = '+str(self.val)+' radians/sec')
            res.append('*DLOAD')
            res.append('%s,CENTRIF,%f,0.,0.,0.,0.,1.,0.' % (cname, radsq))
        elif self.ltype == 'rpm':
            rad = 2.0*pi*(self.val/60.0)
            radsq = rad**2 # square the speed in radians per sec
            res.append('** Rotation force on component: '+cname)
            res.append('** Speed = '+str(self.val)+' rotations/minute')
            res.append('*DLOAD')
            res.append('%s,CENTRIF,%f,0.,0.,0.,0.,1.,0.' % (cname, radsq))
        elif self.ltype == 'matl':
            mat = self.val.name
            res.append('*SOLID SECTION,ELSET='+cname+',MATERIAL='+mat)

        return res


class LinearLoad(object):
    """Makes a load which varies depending on x, y, or z location.

    Load Equation:
        Load(x) = const + mult*(xo - x)

    This load is used to set water pressure, where p = po + rho*g*(xo - x)
    For water pressure this becomes:

        - P(x) = po + (rho*g)*(xo-x)
        - P(x) = po + mult*(xo-x)
        - g is positive, xo > x

    Args:
        ltype (str): load type 'press_fluid' is the only valid option
        comp (Component): component that the load acts on
        const (float): constant term in load equation
        mult (float): mult term in the load equation
        xo (float): xo term in the load equation
        axis (str): axis that the load depends on
            'x' or 'y' are valid options
    """

    def __init__(self, ltype, comp, const, mult, xo, axis='x'):
        self.ltype = ltype
        self.comp = comp
        self.const = const
        self.mult = mult
        self.xo = xo
        self.axis = axis

    def get_val(self, xval):
        """Returns the load function value at the given location.

        Args:
            xval (float): the location where we want to find the load
        """
        res = self.const + (self.xo - xval)*self.mult
        return res

    def get_list(self):
        """Returns a list of [item, val] for plotting."""
        res = []
        faces = self.comp.get_children()
        for face in faces:
            xvals = [getattr(n, self.axis) for n in face.nodes]
            xavg = sum(xvals)/len(xvals)
            pval = self.get_val(xavg)
            res.append([face, pval])
        return res

    def ccx(self):
        """Writes a load for Calculix ccx solver.

        Returns a list strings, where each item is a line to write
        to a Calculix .inp text file.
        """
        res = []
        res.append('** Fluid pressure on component: '+self.comp.name)
        faces = self.comp.get_children()
        for face in faces:
            res.append('*DLOAD')
            xvals = [getattr(n, self.axis) for n in face.nodes]
            xavg = sum(xvals)/len(xvals)
            pval = self.get_val(xavg)
            myline = '%i,P%i,%f' % (face.element.id, face.id, pval)
            res.append(myline)
        return res

class Couple(object):
    """Couples together nodes under items.

    Args:
        comp (Component): component of all the nodes that are coupled
        axes (str): string direction of coupling

            'x': x-motion is identical for all coupled nodes
            'y': y-motion is identical for all coupled nodes
            'xy': x and y motion is identical for all coupled nodes
    """

    def __init__(self, comp, axes='x'):
        self.comp = comp
        self.axes = axes

    def ccx(self):
        """Writes a load for Calculix ccx solver.

        Returns a list strings, where each item is a line to write
        to a Calculix .inp text file.
        """
        res = []
        res.append('** Coupling on component: '+self.comp.name)
        nodes = self.comp.get_children()
        """
        for face in faces:
            res.append('*DLOAD')
            xvals = [getattr(n, self.axis) for n in face.nodes]
            xavg = sum(xvals)/len(xvals)
            pval = self.get_val(xavg)
            myline = '%i,P%i,%f' % (face.element.id, face.id, pval)
            res.append(myline)
        """
        return res
