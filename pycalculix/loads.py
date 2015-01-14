from math import pi

from . import base_classes

class Load(base_classes.Idobj):
    """Makes a load. Many types are supported.
    
    Possible load types:
    fx,fy,fz: forces
    ux,uy,uz: displacements
    thickness: thickness on elements
    matl: matl on elements
    press: pressure
    gravity: gravity
    
    Args:
      ltype (string): load type:
        'fx','fy','fz': force on each axis
        'ux','uy','uz': displacement on each axis
        'thickness': thickness on elements
        'matl': matl on elements
        'press': pressure, + is tension, - is compresion
        'gravity': gravity in x axis direction, - goes towards y axis
        'radps': radians per second rotation
        'rpm': rotations per second rotation
      comp (Component): component that the load is applied to
      val (double or Matl): value of the load
        Matl is needed when setting material loads
        otherwise a double value is used to describe the load
    """
    
    def __init__(self, ltype, comp, val):
        self.ltype = ltype
        self.comp = comp
        self.val = val
        base_classes.Idobj.__init__(self)

    def get_list(self):
        """Returns a list of [item, val] for plotting."""
        res = []
        if self.ltype == 'press':
            # pressure
            faces = self.comp.get_children()
            for f in faces:
                res.append([f,self.val])
        elif self.ltype in ['ux','uy','uz']:
            nodes = self.comp.get_children()
            d = {self.ltype:self.val}            
            for n in nodes:
                res.append( [n,d] )
        return res
        
    def ccx(self):
        """Writes a load for Calculix ccx solver.
        
        Returns a list strings, where each item is a line to write
        to a Calculix .inp text file.        
        """
        res = []
        if self.ltype == 'press':
            # pressure
            res.append('** Pressure on component: '+self.comp.name)
            res.append('*DLOAD')
            faces = self.comp.get_children()
            for f in faces:
                myline = '%i,P%i,%f' % (f.element.id, f.id, self.val)
                res.append(myline)
        elif self.ltype == 'gravity':
            # gravity
            cname = self.comp.name
            res.append('** Gravity on component: '+cname)
            res.append('*DLOAD')
            res.append('%s,GRAV,%f,-1.,0.,0.' % (cname, self.val))            
        elif self.ltype[0] == 'u':
            # displacement
            axis = self.ltype[1]
            anum = {'x':1,'y':2,'z':3}
            axis = anum[axis]
            res.append('*BOUNDARY')
            res.append('%s,%i,%i,%f' % (self.comp.name, axis, axis, self.val))
        elif self.ltype[0] == 'f':
            # force
            res.append('** Force on component: '+self.comp.name)
            nodes = len(self.comp.items)
            fval = self.val/nodes
            axis = self.ltype[1]
            anum = {'x':1,'y':2,'z':3}
            axis = anum[axis]
            res.append('*CLOAD')
            res.append('%s,%i,%f' % (self.comp.name, axis, fval))
        elif self.ltype == 'thickness':
            res.append('*NODAL THICKNESS')
            res.append('%s, %f' % (self.comp.name, self.val))
        elif self.ltype == 'radps':
            radsq = self.val**2 # square the speed in radians per sec
            res.append('** Rotation (radianself, value=radians^2) on component: '+self.comp.name)
            res.append('*DLOAD')
            res.append('%s,CENTRIF,%f,0.,0.,0.,0.,1.,0.' % (self.comp.name, radsq))
        elif self.ltype == 'rpm':
            rad = 2.0*pi*(self.val/60.0)
            radsq = rad**2 # square the speed in radians per sec
            res.append('** Rotation (rpm to radianself, value=radians^2) on component: '+self.comp.name)
            res.append('*DLOAD')
            res.append('%s,CENTRIF,%f,0.,0.,0.,0.,1.,0.' % (self.comp.name, radsq))
        elif self.ltype == 'matl':
            mat = self.val.name
            res.append('*SOLID SECTION,ELSET='+self.comp.name+',MATERIAL='+mat)

        return res


class Load_linear(base_classes.Idobj):
    """Makes a load which varies depending on x, y, or z location.
    
    Load Equation:
      Load(x) = const + mult*(x - xo)

    This load is used to set water pressure, where p = po +rho*g*(x-xo)    
    For water pressure this becomes:
      P(x) = po + (rho*g)*(x-xo)
      P(x) = po + mult*(x-xo)
    
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
        base_classes.Idobj.__init__(self)

    def get_val(self, x):
        """Returns the load function value at the given loaction."""
        res = self.const + (self.xo - x)*self.mult
        return res

    def get_list(self):
        """Returns a list of [item, val] for plotting."""
        res = []
        faces = self.comp.get_children()
        for f in faces:
            xs = [getattr(n, self.axis) for n in f.nodes]
            xavg = sum(xs)/len(xs)
            pval = self.get_val(xavg)
            res.append([f,pval])
        return res

    def ccx(self):
        """Writes a load for Calculix ccx solver.
        
        Returns a list strings, where each item is a line to write
        to a Calculix .inp text file.        
        """
        res = []
        res.append('** Fluid pressure on component: '+self.comp.name)
        faces = self.comp.get_children()
        for f in faces:
            res.append('*DLOAD')
            xs = [getattr(n, self.axis) for n in f.nodes]
            xavg = sum(xs)/len(xs)
            pval = self.get_val(xavg)
            myline = '%i,P%i,%f' % (f.element.id, f.id, pval)
            res.append(myline)    
        return res
