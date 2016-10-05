#include "altitude_hydrostatic.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(altitude_hydrostatic, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<Pressure> p = config_pressure;
  boost::shared_ptr<Temperature> t = config_temperature;
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit height(416, units::m);
  AltitudeHydrostatic a(p, t, lat, height);
  boost::shared_ptr<Altitude> aclone = a.clone();
  blitz::Array<double, 1> expect(19);
  expect =
    47.170180894710057373, 18.825422323132901425, 16.455390829852859724, 11.881822048001820846,
    9.6647632307470789925, 8.1491574419469099411, 7.2074038859224565101, 6.3566118366533359563,
    5.5795506227255309284, 4.8624895301816239979, 4.1960923980757947049, 3.5749978069646335399,
    2.9946606448495125541, 2.4499317590849178927, 1.9364300478578886366, 1.4513850881926142478,
    0.99298497922176331976, 0.55948405200088568989, 0.41600000000000003642;

  for(int i = 0; i < expect.rows(); ++i) {
    BOOST_CHECK_CLOSE(a.altitude(p->pressure_grid()(i)).value.value(), 
                      expect(i),1e-4);
    BOOST_CHECK_CLOSE(aclone->altitude(p->pressure_grid()(i)).value.value(), 
                      expect(i),1e-4);
  }
  expect = 
    9.693291474945075592, 9.7721184143145460865, 9.7793740870791108222, 9.7933884711716387983,
    9.8001974736557020407, 9.8048535643795133865, 9.8077476452213119273, 9.8103629665601452814,
    9.8127525378688851276, 9.8149583482049429506, 9.8170089478969906338, 9.8189207203840105365,
    9.8207075494057338716, 9.8223851818337450936, 9.8239670247695354277, 9.8254615511560636776,
    9.8268742924496201852, 9.8282106958417934095, 9.8286531264705985222;
  for(int i = 0; i < expect.rows(); ++i) {
    BOOST_CHECK_CLOSE(a.gravity(p->pressure_grid()(i)).value.value(), 
                      expect(i),1e-4);
    BOOST_CHECK_CLOSE(aclone->gravity(p->pressure_grid()(i)).value.value(), 
                      expect(i),1e-4);
  }

  // Turn on for debugging
  if (false) {
    std::cerr << setprecision(20);
    std::cerr << "Altitude: " << expect.rows() << std::endl;
    for(int i = 0; i < expect.rows(); ++i) {
      std::cerr << a.altitude(p->pressure_grid()(i)).value.value() << ", ";
    }
    std::cerr << std::endl;
    std::cerr << "Gravity: " << expect.rows() << std::endl;
    for(int i = 0; i < expect.rows(); ++i) {
      std::cerr << a.gravity(p->pressure_grid()(i)).value.value() << ", ";
    }
    std::cerr << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(jacobian)
{
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit height(416, units::m);
  AltitudeHydrostatic a(config_pressure, config_temperature, 
                        lat, height);
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());

  ArrayAdWithUnit<double, 1> p(config_pressure->pressure_grid());
  Array<AutoDerivative<double>, 1> alt(p.rows());
  Array<AutoDerivative<double>, 1> grav(p.rows());
  Array<double, 1> alt0(alt.shape()), grav0(grav.shape());
  for(int i = 0; i < p.rows(); ++i) {
    alt(i) = a.altitude(p(i)).value;
    grav(i) = a.gravity(p(i)).value;
    alt0(i) = alt(i).value();
    grav0(i) = grav(i).value();
  }
  Array<double, 2> altjac = FullPhysics::jacobian(alt);
  Array<double, 2> gravjac = FullPhysics::jacobian(grav);
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 1> altjacfd(alt0.shape()), gravjacfd(grav0.shape());
    for(int j = 0; j < p.rows(); ++j) {
      altjacfd(j) = (a.altitude(p(j)).value.value() - alt0(j)) / epsilon(i);
      gravjacfd(j) = (a.gravity(p(j)).value.value() - grav0(j)) / epsilon(i);
    }
    if(false) {                  // Can turn this off to dump values,
                                // if needed for debugging
      double diff = max(abs(altjac(Range::all(), i) - altjacfd));
      std::cerr.precision(20);
      if(diff > 0)
        std::cerr << i << ": " << altjac(Range::all(), i) << "\n"
                  << altjacfd << "\n"
                  << diff << "\n";
      diff = max(abs(gravjac(Range::all(), i) - gravjacfd));
      if(diff > 0)
        std::cerr << i << ": " << gravjac(Range::all(), i) << "\n"
                  << gravjacfd << "\n"
                  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(altjac(Range::all(), i), altjacfd, 1e-8);
    BOOST_CHECK_MATRIX_CLOSE_TOL(gravjac(Range::all(), i), gravjacfd, 1e-8);
  }
}

BOOST_AUTO_TEST_SUITE_END()
