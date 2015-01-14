from . import base_classes

class Material(base_classes.Idobj):
    """Makes a linear elastic meterial.

    Args:
      name (str): a unique matl name

    Attributes:
      id (int): material id number
      density (float): density in mass/volume units
      pratio (float): poisson's ratio (unitless)
      youngs (float): young's modulus in Force/area = stress units
      conductivity (float): Power / (distance-temp)
      spec_heat (float): specific heat (energy/(mass-temp)
      thermal_exp (dict): a dict storing temperature dependent thermal
        expansion properties
        dict['data'] = zip of (temp, thermal_expansion)
          Thermal expansion is in strain per temperature
        dict['tzero'] = the temperature zero point
    """

    def __init__(self, name):
        self.name = name
        base_classes.Idobj.__init__(self)        

    def set_mech_props(self, density, youngs, pratio):
        """Sets the mechanical properties: density, youngs, poisson_ratio.

        Args:
          density (float): density in mass/volume units
          pratio (float): poisson's ratio (unitless)
          youngs (float): young's modulus in Force/area = stress units
        """
        self.density = density
        self.youngs = youngs
        self.pratio = pratio

    def set_therm_props(self, conductivity, spec_heat):
        """Sets the thermal properties: conductivity, specifc_heat.

        Args:
          conductivity (float): Power / (distance-temp)
          spec_heat (float): specific heat (energy/(mass-temp)
        """
        self.conductivity = conductivity
        self.spec_heat = spec_heat

    def set_therm_expan(self, temps, alphas, tzero):
        """Sets the thermal expansion of the material.

        Args:
          temps (list): list of temperatures
          alphas (list): list of thermal expansion alphas
                         length must be the same as the passed temps
          tzero (float): temperature zero point
        """
        self.thermal_exp = {}
        self.thermal_exp['data'] = zip(alphas, temps)
        self.thermal_exp['tzero'] = tzero

    def ccx(self):
        """Returns a list of text strings for ccx defining the material."""
        res = []
        res.append('*MATERIAL,NAME='+self.name)
        if hasattr(self,'youngs'):
            res.append('*ELASTIC')
            res.append(str(self.youngs)+','+str(self.pratio))
        if hasattr(self, 'density'):
            res.append('*DENSITY')
            res.append(str(self.density))
        if hasattr(self, 'conductivity'):
            res.append('*CONDUCTIVITY')
            res.append(str(self.conductivity))
        if hasattr(self, 'spec_heat'):
            res.append('*SPECIFIC HEAT')
            res.append(str(self.spec_heat))
        if hasattr(self, 'thermal_exp'):
            res.append('*EXPANSION,ZERO='+str(self.tzero))
            for pair in self.thermal_exp:
                res.append(str(pair[0])+' '+str(pair[1]))
        return res