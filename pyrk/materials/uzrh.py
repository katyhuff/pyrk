from utilities.ur import units
from density_model import DensityModel
from materials.material import Material


class UZrH(Material):
    """This class represents Uranium-Zirconium hydride nuclear reactor fuel
    ($$U{0.3}ZrH_{1.6}$$). It inherits from the material class and possesses
    attributes intrinsic to $$U{0.3}ZrH_{1.6}$$
    Note that this notation $$U_yZrH_x$$ indicates that x is the H/Zr atom
    ratio and y is the U/Zr atom ratio.

    References:
        D. Olander, E. Greenspan, H. D. Garkisch, and B. Petrovic,
        “Uranium–zirconium hydride fuel properties,” Nuclear Engineering and
        Design, vol. 239, no. 8, pp. 1406–1424, Aug. 2009.
    """
    def __init__(self,
                 name="uzrh",
                 x=1.6,
                 e=0.20,
                 w_u=0.45*units.g/(units.cm, 3)):
        """Initalizes a material

        :param name: The name of the component (i.e., "fuel" or "cool")
        :type name: str.
        """
        Material.__init__(self,
                          name=name,
                          k=self.thermal_conductivity(),
                          cp=self.specific_heat_capacity(),
                          dm=self.density())
        self.x = x
        self.e = e
        self.w_u = w_u
        self.mu = self.m_u(self.e)
        self.y = self.y(self.mu, w_u)

    def thermal_conductivity(self):
        """UZrH thermal conductivity in [W/m-K]
        """
        return 18.0*units.watt/(units.meter*units.kelvin)

    def specific_heat_capacity(self):
        """Specific heat capacity of uzrh [J/kg/K]
        TODO: varies with temperature in J/mole-K :
        $$(25 + 4.7x) + (0.31+2.01x)T/100 + (1.9+6.4x)/T^2 × 10^{−5}$$
        """
        return 4181.3*units.joule/(units.kg*units.kelvin)

    def m_u(self,
            e=0.20):
        """M_u is the atomit wieght of the uranium mixture

        :param e: is the enrichment of u235
        :type e: float (3% = 0.03)
        """
        return 235*e + 238*(1-e)

    def y(self, mu, w_u):
        """U/Zr atom ratio

        :param mu: M_u is the atomic weight of the uranium mix
        :type mu: float
        :param w_u: is the weight fraction of uranium
        :type w_u: float
        """
        y = (91.2/mu)*(w_u/(1-w_u))
        return y

    def rho_ZrH(self,
                x):
        """$$\rho_{ZrH}$$ The density of the ZrH in the fuel
        :param x: is the atom ratio of H to Zr
        :type x: float
        """
        if x < 1.6:
            rho = ((0.154 + 0.0145*x) - 1)*units.g/pow(units.cm, 3)
        elif x >= 1.6:
            rho = ((0.171 + 0.0042*x) - 1)*units.g/pow(units.cm, 3)

        to_ret = rho.to(units.kg/pow(units.m, 3))
        return to_ret

    def rho_UZrH(self,
                 w_u):
        """
        w_u

        """
        rho_U0 = 19.9*units.g/pow(units.cm, 3)
        denom = w_u/rho_U0 + (1-w_u)/self.rho_ZrH(self.x)
        rho = 1/denom
        return rho.to(units.kg/pow(units.m, 3))

    def density(self):
        """UZrH density as a funciton of T. [kg/m^3]
        TODO: Pressure-based density function (tait) for uzrh.

        """
        return DensityModel(a=self.rho_UZrH(self.w_u),
                            model="constant")
