from nose.tools import *
from full_physics import *
import os
import numpy as np
import numpy.testing as nptest
import datetime;
import time;

#def test_solar_doppler_shift_polynomial():
# Working on solar stuff, so for now don't run this test.
if(False):
    # Time as a datetime. Note you can also use a float as a unix time
    # stamp if that is more convenient. We include fractions of a second

    tm = datetime.datetime(2006,9,14,12,27,22,100)

    # Somewhat oddly, python doesn't support timezones out of the box.
    # Rather than depending on a third party library for a unit test,
    # we just convert from local time to GMT time

    tm = tm - datetime.timedelta(seconds=time.timezone)
    if(time.daylight):
        tm += datetime.timedelta(hours=1)

    p = SolarDopplerShiftPolynomial(tm,
                                    DoubleWithUnit(77.1828918457,Unit("deg")),
                                    DoubleWithUnit(74.128288269,Unit("deg")),
                                    DoubleWithUnit(167.495071411,Unit("deg")),
                                    DoubleWithUnit(416,Unit("m")))

    assert_almost_equal(p.solar_distance.value, 1.0060305651331354)

    wn = SpectralDomain([12929.94, 12979.93, 13029.93, 13079.93, 
                         13129.93, 13179.93])
    expected = np.array([12929.919173650407,
                         12979.909093131146,
                         13029.909012595777,
                         13079.908932060409,
                         13129.908851525041,
                         13179.908770989672])
    nptest.assert_array_almost_equal(p.doppler_stretch(wn).data, expected)
    
