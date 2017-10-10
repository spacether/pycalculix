"""This module stores material classes.
"""

from . import base_classes

class Material(base_classes.Idobj):
    """Makes a linear elastic meterial.

    Args:
        name (str): a unique matl name

    Attributes:
        - id (int): material id number
        - mechtype (str): defines mechanical material type: linear or non-linear material
        - density (float): density in mass/volume units
        - pratio (float): poisson's ratio (unitless)
        - youngs (float): young's modulus in Force/area = stress units
        - mechtpye (str): mechanical material type
            options: 'linear' or 'nonlinear'
        - exponent (float): exponent of Ramberg-Osgood stress-strain equation
        - yield_stress (float): yield stress of material
        - yield_offset (float): yield offset of Ramberg-Osgood stress-strain equation
        - conductivity (float): thermal conductivity, Power / (distance-temp)
        - spec_heat (float): specific heat (energy/(mass-temp)
        - thermal_exp (dict): a dict storing temperature dependent thermal
            -- expansion properties

            -- Thermal expansion is in strain per temperature

            - dict['data'] = zip of (temp, thermal_expansion)
            - dict['tzero'] = the temperature zero point
    """

    def __init__(self, name):
        self.name = name
        base_classes.Idobj.__init__(self)
        # mechanical
        self.mechtype = None
        self.density = None
        self.youngs = None
        self.pratio = None
        self.exponent = None
        self.yield_stress = None
        self.yield_offset = None
        # thermal
        self.conductivity = None
        self.spec_heat = None
        # thermal growth
        self.thermal_exp = {}

    def set_mech_props(self, density, youngs, pratio, mechtype='linear', exponent=0., yield_stress=0., yield_offset=0.):
        """Sets the mechanical properties: density, youngs, poisson_ratio.

        Args:
          - density (float): density in mass/volume units
          - pratio (float): poisson's ratio (unitless)
          - youngs (float): young's modulus in Force/area = stress units
          - mechtpye (str): mechanical material type
            options: 'linear' or 'nonlinear'
          - exponent (float): exponent of Ramberg-Osgood stress-strain equation
          - yield_stress (float): yield stress of material
          - yield_offset (float): yield offset of Ramberg-Osgood stress-strain equation
        """
        self.density = density
        self.youngs = youngs
        self.pratio = pratio
        self.mechtype = mechtype
        self.exponent = exponent
        self.yield_stress = yield_stress
        self.yield_offset = yield_offset

    def set_therm_props(self, conductivity, spec_heat):
        """Sets the thermal properties: conductivity, specifc_heat.

        Args:
          - conductivity (float): Power / (distance-temp)
          - spec_heat (float): specific heat (energy/(mass-temp)
        """
        self.conductivity = conductivity
        self.spec_heat = spec_heat

    def set_therm_expan(self, alphas, temps=None, tzero=None):
        """Sets the thermal expansion of the material.

        Args:
          - alphas (float or list): list of thermal expansion alphas
                         length must be the same as the passed temps
          - temps (list): list of temperatures
          - tzero (float): temperature zero point
        """
        self.thermal_exp = {}
        isnum = (isinstance(alphas, float) or isinstance(alphas, int))
        if isinstance(alphas, list) and isinstance(temps, list):
            self.thermal_exp['data'] = zip(temps, alphas)
        elif temps == None and isnum:
            self.thermal_exp['data'] = [alphas]
        if tzero != None:
            self.thermal_exp['tzero'] = tzero

    def ccx(self):
        """Returns a list of text strings for ccx defining the material."""
        res = []
        res.append('*MATERIAL,NAME='+self.name)

        if self.mechtype == 'linear':
            res.append('*ELASTIC')
            res.append(str(self.youngs)+','+str(self.pratio))
        elif self.mechtype == 'nonlinear':
            res.append('*DEFORMATION PLASTICITY')
            res.append(','.join([str(self.youngs),str(self.pratio),str(self.yield_stress),str(self.exponent),str(self.yield_offset)]))
        if self.density != None:
            res.append('*DENSITY')
            res.append(str(self.density))
        if self.conductivity != None:
            res.append('*CONDUCTIVITY')
            res.append(str(self.conductivity))
        if self.spec_heat != None:
            res.append('*SPECIFIC HEAT')
            res.append(str(self.spec_heat))
        if self.thermal_exp != {}:
            if 'tzero' in self.thermal_exp:
                res.append('*EXPANSION,ZERO='+str(self.thermal_exp['tzero']))
            else:
                res.append('*EXPANSION')
            for pair in self.thermal_exp['data']:
                if len(pair) == 1:
                    res.append(str(pair[0]))
                else:
                    # temp, val
                    res.append('%f %f' % (pair[0], pair[1]))
        return res
