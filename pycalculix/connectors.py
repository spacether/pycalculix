"""This module stores connector classes, like Contact.
"""
from . import base_classes

class SurfaceInteraction(base_classes.Idobj):
    """Makes a surface interaction object.

    Args:
        int_type (str): interaction type

            * 'EXPONENTIAL'
            * 'LINEAR'

        *args: following arguments

            * int_type = 'EXPONENTIAL'
                * c0, then p0 must be passed

            * int_type = 'LINEAR'
                * k must be passed

    """

    def __init__(self, int_type, *args):
        self.int_type = int_type
        self.c0 = None
        self.p0 = None
        self.k = None
        if self.int_type == 'EXPONENTIAL':
            self.c0 = args[0]
            self.p0 = args[1]
        elif self.int_type == 'LINEAR':
            self.k = args[0]
        base_classes.Idobj.__init__(self)

    @property
    def name(self):
        """SurfaceInteraction name."""
        return 'SI%i' % self.id

    def ccx(self):
        """Writes the surface interaction for the ccx solver."""
        res = ["*SURFACE INTERACTION,NAME=%s" % self.name]
        res.append("*SURFACE BEHAVIOR,PRESSURE-OVERCLOSURE=%s" % self.int_type)
        if self.int_type == 'EXPONENTIAL':
            res.append("%f,%f" % (self.c0, self.p0))
        elif self.int_type == 'LINEAR':
            res.append("%e" % self.k)
            #res.append("%e,1." % self.k)
        return res

class Contact(base_classes.Idobj):
    """Makes a contact which will be between lines which have faces on them.

    Args:
        master_comp (Component): component of master element faces
        slave_comp (Component): component of slave element faces
        surf_int (SurfaceInteraction): object which stores closure behavior
        surf_to_surf (bool): if True surface to surface is used, if False node
            to surface is used
    """

    def __init__(self, master_comp, slave_comp, surf_int, surf_to_surf=True):
        self.master_comp = master_comp
        self.slave_comp = slave_comp
        self.surf_int = surf_int
        self.surf_to_surf = surf_to_surf
        base_classes.Idobj.__init__(self)

    @property
    def name(self):
        """Contact name."""
        return 'CONT%i' % self.id

    def ccx(self):
        """Writes the contact pair for the ccx solver."""
        res = []
        end = ''
        int_name = self.surf_int.name
        m_name = self.master_comp.name
        s_name = self.slave_comp.name
        if self.surf_to_surf:
            end = ',TYPE=SURFACE TO SURFACE'
        line = "*CONTACT PAIR,INTERACTION=%s%s" % (int_name, end)
        res.append(line)
        line = "%s,%s" % (s_name, m_name)
        res.append(line)
        return res
