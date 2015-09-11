from utilities.ur import units


class HeatRemoval(object):
    """This is the default heat removal class from which all others are
    derived.  The default is constant, 100% efficiency such that the inlet
    temperature is constant throughout the simulation..


    | T = T_in     __________________________________
    |
    |
    |            t0                                tf
    """

    def __init__(self, timer,
                 T_in=0.0*units.kelvin,
                 delay=0.0*units.seconds):
        """Creates a HeatRemoval object to represent a heat sink.

        :param T_in: The design point inlet temperature (100% efficiency)
        :type T_in: float (units of kelvin)
        :param timer: The timer object for the simulation
        :type timer: Timer
        """
        self.T_in = T_in.to('kelvin')
        self.timer = timer
        self.delay = delay.to('seconds')
        self.vals = [self.f(t_idx) for t_idx in range(timer.timesteps())]

    def f(self, x):
        """ A function that returns the cooled inlet temperature at each timestep.

        :param x: The timestep index
        :type x: int
        """
        return self.T_in.to('kelvin')

    def T(self, t_idx, T_out):
        """
        Return the temperature T
        """
        return self.T_in.to('kelvin')


class StepHeatRemoval(HeatRemoval):
    """Returns a Heaviside step function.
    The default for this model is a complete loss of heat removal.


    | eff_init _________________
    |                           |
    |                           |
    |                           |
    |                           |
    |                           |
    |                           |
    | eff_final                 |_____________________
    |
    |                       t_step
    """
    def __init__(self,
                 timer,
                 T_in=1.0*units.kelvin,
                 t_step=1.0*units.seconds,
                 eff_init=1.0,
                 eff_final=0.0):
        """Returns a Heaviside step function.
        Default is a 100% loss of heat sink at time t=1s.

        :param timer: The timer object for the simulation
        :type timer: Timer
        :param T_in: The design point inlet temperature (100% efficiency)
        :type T_in: float (units of kelvin)
        :param t_step: The time at which the step (or drop) occurs
        :type t_step: float (units seconds)
        :param eff_init: Initial efficiency of cooling
        :type eff_init: float
        :param eff_final: Final efficiency of cooling
        :type eff_final: float
        """
        self.eff_init = eff_init
        self.eff_final = eff_final
        self.t_step = t_step.to('seconds')
        HeatRemoval.__init__(self, timer=timer)

    def f(self, x):
        """ A function that returns the cooled inlet temperature at each timestep.

        :param x: The timestep index
        :type x: int
        """
        if x < self.timer.t_idx(self.t_step):
            return self.eff_init
        else:
            return self.eff_final

    def T(self, t_idx, T_out):
        """
        Return the temperature T
        """
        return T_out*(1 - self.vals[t_idx]) + self.T_in*self.vals[t_idx]
