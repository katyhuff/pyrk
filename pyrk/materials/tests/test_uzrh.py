from nose.tools import assert_equal, assert_almost_equal, assert_true, \
    assert_false, assert_raises, assert_is_instance, with_setup

from materials import uzrh
from utilities.ur import units

name = "testname"
tester = uzrh.UZrH(name=name)


T0 = 500.0*units.kelvin
k_uzrh = 18.0*units.watt/(units.meter*units.kelvin)
cp_uzrh = 4181.3*units.joule/(units.kg*units.kelvin)
rho_at_temp_zero = 2413.2172*units.kg/units.meter**3
w_u = 0.45


def test_constructor():
    assert_equal(tester.name, name)
    assert_equal(tester.k, k_uzrh)
    assert_equal(tester.cp, cp_uzrh)
    rho_at_temp_zero = tester.rho(T0)
    # currently, rho is constant. when that changes, change test
    assert_equal(tester.rho(T0+100*units.kelvin), rho_at_temp_zero)
    assert_equal(tester.w_u, w_u)
