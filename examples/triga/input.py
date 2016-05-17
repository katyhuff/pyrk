from utilities.ur import units
import th_component as th
import math
from materials.uzrh import UZrH
from materials.water import Water
from timer import Timer
from reactivity_insertion import ImpulseReactivityInsertion

#############################################
#
# User Workspace
#
#############################################

# Thermal hydraulic params
# Temperature feedback of reactivity
# A. L. Costa, C. Pereira, P. A. L. Reis, M. A. F. Veloso, A. Z. Mesquita,
# and H. V. Soares, "Thermal hydraulic analysis of the IPR-R1 TRIGA research
# reactor using a RELAP5 model," Nuclear Engineering and Design, vol. 240,
# no. 6, pp. 1487-1494, Jun. 2010.
alpha_f = (-1.1*units.pcm/units.kelvin)
# below from random guess
t_fuel = 500.00*units.kelvin
t_cool = 500.00*units.kelvin

kappa = 0.00  # TODO if you fix omegas, kappa ~ 0.06

# Initial time
t0 = 0.00*units.seconds

# Timestep
dt = 0.0005*units.seconds

# Final Time
tf = 10.0*units.seconds

# Time at which temperature feedbacks start
t_feedback = 0.5*units.seconds


def vol_cyl(r, h):
    return (math.pi*pow(r, 2)*h).to(units.meter**3)

# volumes
# P. A. L. Reis, A. L. Costa, C. Pereira, M. A. F. Veloso, A. Z. Mesquita, H. V.
# Soares, and G. de P. Barros, "Assessment of a RELAP5 model for the IPR-R1
# TRIGA research reactor," Annals of Nuclear Energy, vol. 37, no. 10, pp.
# 1341-1350, Oct. 2010.
#
# A. L. Costa, P. A. Reis, C. A. Silva, C. Pereira, M. A. F. Veloso, B. T.
# Guerra, H. V. Soares, and A. Z. Mesquita, "Safety studies and general
# simulations of research reactors using nuclear codes," Nuclear Power-System
# Simulations and Operation, 2011.
#
# A. Z. Mesquita and H. C. Rezende, "Experimental heat transfer analysis of the
# IPR-R1 TRIGA reactor," in Proceedings of the 3rd World TRIGA USERS Conference.
# Centro De Desenvolvimento Da Tecnologia Nuclear-CDTN/CNEN, Brazil, 2006.

d_rod = 35.6*units.millimeter
r_rod = d_rod/2
h_rod = 355.6*units.millimeter
n_fuel_rods = 63
v_fuel_rod = vol_cyl(r_rod, h_rod)
v_fuel = n_fuel_rods*v_fuel_rod
r_pool = (1.92/2.)*units.meters  # http://cdn.intechopen.com/pdfs-wm/18106.pdf
h_pool = 6.625*units.meters  # http://cdn.intechopen.com/pdfs-wm/18106.pdf
v_pool = vol_cyl(r_pool, h_pool)
v_water = v_pool - v_fuel

h_water = 7.0*units.kilowatt/pow(units.meter, 2)/units.kelvin
a_water = 2.0*math.pi*r_rod*h_rod*n_fuel_rods  # tot surface area

#############################################
#
# Required Input
#
#############################################

# Total power, Watts, thermal
# Here we assume 1 fuel pin for simplicity
power_tot = 250.0*units.kilowatt

# Timer instance, based on t0, tf, dt
ti = Timer(t0=t0, tf=tf, dt=dt, t_feedback=t_feedback)

# Number of precursor groups
n_pg = 6

# Number of decay heat groups
n_dg = 0

# Fissioning Isotope
fission_iso = "u235"

# Spectrum
spectrum = "thermal"

# Feedbacks, False to turn reactivity feedback off. True otherwise.
feedback = True

# External Reactivity
rho_ext = ImpulseReactivityInsertion(timer=ti,
                                     t_start=1.0*units.seconds,
                                     t_end=2.0*units.seconds,
                                     rho_init=0.0*units.delta_k,
                                     rho_max=0.001*units.delta_k)

# maximum number of internal steps that the ode solver will take
nsteps = 1000000

fuel = th.THComponent(name="fuel",
                      mat=UZrH(name="fuel"),
                      vol=v_fuel,
                      T0=t_fuel,
                      alpha_temp=alpha_f,
                      timer=ti,
                      heatgen=True,
                      power_tot=power_tot)

cool = th.THComponent(name="cool",
                      mat=Water(name="water"),
                      vol=v_water,
                      T0=t_cool,
                      timer=ti)

# The clad convects with the coolant
fuel.add_convection('cool', h=h_water, area=a_water)
cool.add_convection('fuel', h=h_water, area=a_water)

components = [fuel, cool]
