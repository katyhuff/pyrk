from utilities.ur import units
from density_model import DensityModel
from materials.material import Material


class Water(Material):
    """This class represents Water. It inherits from the material
    class and possesses attributes intrinsic to water.

    References:
    https://en.wikipedia.org/wiki/Properties_of_water
    "Thermal Conductivity of some common Materials". The Engineering ToolBox. Retrieved 2011-11-22.
    """
    def __init__(self, name="water"):
        """Initalizes a material

        :param name: The name of the component (i.e., "fuel" or "cool")
        :type name: str.
        """
        Material.__init__(self,
                          name=name,
                          k=self.thermal_conductivity(),
                          cp=self.specific_heat_capacity(),
                          dm=self.density())

    def thermal_conductivity(self):
        """Water thermal conductivity in [W/m-K]
        """
        return 0.58*units.watt/(units.meter*units.kelvin)

    def specific_heat_capacity(self):
        """Specific heat capacity of water [J/kg/K]
        """
        return 4181.3*units.joule/(units.kg*units.kelvin)

    def density(self):
        """
        Water density as a funciton of T. [kg/m^3]
        TODO: Pressure-based density function (tait) for water.

        """
        return DensityModel(a=1000*units.kg/(units.meter**3),
                            model="constant")

