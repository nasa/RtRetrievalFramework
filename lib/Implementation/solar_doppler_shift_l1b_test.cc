#include "solar_doppler_shift_l1b.h"
#include "unit_test_support.h"
#include "fp_exception.h"
#include "old_constant.h"
#include "default_constant.h"
#include "fts_run_log.h"

using namespace FullPhysics;
using namespace boost::posix_time;
using namespace boost::gregorian;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(solar_doppler_shift_l1b, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  DoubleWithUnit solar_distance(1.0060305651331354, OldConstant::AU);
  DoubleWithUnit solar_velocity(-0.000278885, 
				OldConstant::AU / units::day);
  SolarDopplerShiftL1b p(solar_distance, solar_velocity);
  BOOST_CHECK_CLOSE(p.solar_distance().convert(OldConstant::AU).value, 
                    1.0060305651331354, 1e-3); 
  Array<double, 1> wn(6);
  wn = 12929.940000000001,
    12979.930000000000,
    13029.930000000000,
    13079.930000000000,
    13129.930000000000,
    13179.930000000000;
  Array<double, 1> res = p.doppler_stretch(wn).wavenumber();
  BOOST_CHECK_CLOSE(res(0), 12929.919173650407, 1e-3);
  BOOST_CHECK_CLOSE(res(1), 12979.909093131146, 1e-3);
  BOOST_CHECK_CLOSE(res(2), 13029.909012595777, 1e-3);     
  BOOST_CHECK_CLOSE(res(3), 13079.908932060409, 1e-3);     
  BOOST_CHECK_CLOSE(res(4), 13129.908851525041, 1e-3);     
  BOOST_CHECK_CLOSE(res(5), 13179.908770989672, 1e-3);     
}

BOOST_AUTO_TEST_SUITE_END()
