#include "rayleigh.h"
#include "pressure.h"
#include "configuration_fixture.h"
#include "altitude_hydrostatic.h"
#include "atmosphere_oco.h"
#include "unit_test_support.h"
#include "hdf_file.h"
#include "default_constant.h"
using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(rayleigh, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(cross_section)
{
  DoubleWithUnit c = Rayleigh::cross_section(DoubleWithUnit(532, units::nm));
  BOOST_CHECK_CLOSE(c.convert(Unit("m^2")).value, 5.14127e-31, 1e-4);
}

BOOST_AUTO_TEST_CASE(basic)
{
  // Will replace this
  // HdfFile h(test_data_dir() + "l2_fixed_level_static_input.h5");
  // double depolar_fact = h.read_attribute<double>("/Gas/Air/depolarization_factor");
  // double wdc_a = h.read_attribute<double>("/Earth/wavelength_dependence_coefficient_a");
  // double wdc_b =
  // h.read_attribute<double>("/Earth/wavelength_dependence_coefficient_b");
  DefaultConstant constant;
  boost::shared_ptr<Pressure> p = config_pressure;
  boost::shared_ptr<Temperature> t = config_temperature;
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit height(416, units::m);
  std::vector<boost::shared_ptr<Altitude> > alt;
  alt.push_back(boost::shared_ptr<Altitude>(new AltitudeHydrostatic(p, t, 
				    lat, height)));
  Rayleigh r(p, alt, constant);
  Array<double, 1> od_expect(18);
  od_expect = 0.0016576665656406141224, 0.00071940274937751041116,
    0.002395191514237895395, 0.0019142978250297388882,
    0.0016740605906486161723, 0.0011953048701276952712,
    0.0011949700946174662657, 0.0011946660498398202488,
    0.0011943869014422443895, 0.0011941282697030431315,
    0.0011938874856637184726, 0.001193662889602569192,
    0.0011934526013696337889, 0.0011932547696655935576, 
    0.0011930680393243638972, 0.0011928916227435914631,
    0.0011927248631764474605, 0.00040945552379681356961;
  BOOST_CHECK_MATRIX_CLOSE_TOL(r.optical_depth_each_layer(12929.94, 0).value(), 
			       od_expect, 1e-6);
  od_expect = 0.0016578529204295140969, 0.0007194836245970951099,
    0.0023954607815430718765, 0.0019145130302914421631,
    0.0016742487884529530257, 0.0011954392462388873167,
    0.0011951044330932129867, 0.0011948003541348687674,
    0.0011945211743554446163, 0.001194262513540876678,
    0.0011940217024326228897, 0.0011937970811223960407,
    0.0011935867692488685496, 0.0011933889153045984782,
    0.0011932021639711536901, 0.0011930257275576383862,
    0.0011928589492433926394, 0.00040950155476545508585;
  BOOST_CHECK_MATRIX_CLOSE_TOL(r.optical_depth_each_layer(12930.30, 0).value(), 
			       od_expect, 1e-6);
}

BOOST_AUTO_TEST_CASE(jacobian)
{
  Rayleigh& r = *(dynamic_cast<const AtmosphereOco&>(*config_atmosphere).rayleigh_ptr());
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());

  ArrayAd<double, 1> rgrid = r.optical_depth_each_layer(12929.94, 0);
  Array<double, 1> rgrid0(rgrid.shape());
  rgrid0 = rgrid.value();
  Array<double, 2> jac = rgrid.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 1> jacfd(rgrid0.shape());
    jacfd = (r.optical_depth_each_layer(12929.94, 0).value() - rgrid0)
      / epsilon(i);
    if(false) {			// Can turn this off to dump values,
				// if needed for debugging
      double diff = max(abs(jac(Range::all(), i) - jacfd));
      if(diff > 0)
	std::cerr << i << ": " << jac(Range::all(), i) << "\n"
		  << jacfd << "\n"
		  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), i), jacfd, 1e-11);
  }
}

BOOST_AUTO_TEST_SUITE_END()
