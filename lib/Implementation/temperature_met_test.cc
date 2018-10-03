#include "temperature_met.h"
#include "atmosphere_fixture.h"
#include "acos_sounding_id.h"
#include "unit_test_support.h"
#include "acos_met_file.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(temperature_met, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::string sid = "20091009203401";
  HdfFile sfile(test_data_dir() + "in/sounding_id.h5");
  std::vector<boost::shared_ptr<HdfSoundingId> > sidv = 
    AcosSoundingId::create(sfile, sid);
  boost::shared_ptr<AcosMetFile> e(new AcosMetFile(test_data_dir() + "in/ecmwf.h5", 
					       sidv[0], sidv.size() > 1));
  boost::shared_ptr<Pressure> p = config_pressure;
  StateVector sv;
  TemperatureMet t1(e, p, 0, true);
  sv.add_observer(t1);
  Array<double, 1> temp_expect(19);
  temp_expect = 244.2, 214.553, 218.029, 222.544, 218.341, 221.37, 227.38,
    233.493, 239.376, 244.52, 248.708, 251.979, 254.537, 256.655, 258.521,
    260.155, 261.747, 261.732, 258.598;
  for(int i = 0; i < temp_expect.rows(); ++i)
    BOOST_CHECK_CLOSE(t1.temperature(p->pressure_grid()(i)).convert(units::K).value.value(),
     		      temp_expect(i), 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

// BOOST_FIXTURE_TEST_SUITE(temperature_fixed_level_jac, ConfigurationFixture)

// BOOST_AUTO_TEST_CASE(jacobian)
// {
//   Temperature& t = *config_temperature;
//   StateVector& sv = *config_state_vector;
//   Array<double, 1> sv0(sv.state().copy());

//   // The bottom two levels of the temperature happen to be
//   // identical. That means that when we are between these two levels,
//   // we get no change in temperature as the pressure changes. To make
//   // sure we test this dependency, move the surface pressure so we are
//   // at a higher level where there is a temperature dependency. We use
//   // the fact that the surface pressure is index 21 in the state vector.
//   const int surface_pressure_index = 21;
//   sv0(surface_pressure_index) = 83000;
//   sv.update_state(sv0);
  
//   ArrayAd<double, 1> tgrid = t.temperature_grid();
//   Array<double, 1> tgrid0(tgrid.rows());
//   tgrid0 = tgrid.value();
//   Array<double, 2> jac(tgrid.jacobian().copy());
//   for(int i = 0; i < sv.state().rows(); ++i) {
//     Array<double, 1> svn(sv0.copy());
//     svn(i) += epsilon(i);
//     sv.update_state(svn);
//     Array<double, 1> jacfd(tgrid0.shape());
//     jacfd = (t.temperature_grid().value() - tgrid0) / epsilon(i);
//     if(false) {			// Can turn this off to dump values,
// 				// if needed for debugging
//       std::cerr.precision(20);
//       double diff = max(abs(jac(Range::all(), i) - jacfd));
//       if(diff > 0)
// 	std::cerr << i << ": " << jac(Range::all(), i) << "\n"
// 		  << jacfd << "\n"
// 		  << diff << "\n";
//     }
//     BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), i), jacfd, 1e-3);
//   }
// }

// BOOST_AUTO_TEST_SUITE_END()

// BOOST_FIXTURE_TEST_SUITE(temperature_fixed_level_ecmwf, ConfigurationEcmwfFixture)

// BOOST_AUTO_TEST_CASE(ecmwf)
// {
//   IfstreamCs expected_data(test_data_dir() + "expected/absorber/ecmwf");
//   const TemperatureFixedLevel& tecmwf = 
//     dynamic_cast<TemperatureFixedLevel&>(*config_temperature);
//   Array<double, 1> texpect(20);
//   texpect =
//     243.90097378304642461, 214.55495206293338128, 218.03548301856420721,
//     222.54439025512246531, 218.33300494385892421, 221.37919630741211563,
//     227.3962815102570687,
//     233.51051500723281151, 239.38876092919977623, 244.52113267396558172,
//     248.71471551251292453, 251.98636348202509794, 254.54141229243495559,
//     256.65775880671048981, 
//     258.52334065260964735, 260.15648388783131395, 261.74845655310156189,
//     261.7317782935307946, 256.31781429833779384, 257.36699412179075352;
//   BOOST_CHECK_MATRIX_CLOSE(tecmwf.temperature_levels().value(), texpect);
// }

// BOOST_AUTO_TEST_SUITE_END()
