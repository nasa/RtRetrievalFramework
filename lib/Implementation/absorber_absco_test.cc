#include "absorber_absco.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"
#include "atmosphere_fixture.h"
#include "lua_state.h"
#include <cstdlib>

using namespace FullPhysics;
using namespace blitz;

// Fixture that uses our tccon sounding 1 test. Right now this is only
// used here, we can move this out to its own file if we end up using
// this in other places.
class TcconSounding1Fixture: public GlobalFixture {
private:
  static boost::shared_ptr<LuaState> ls;
  blitz::Array<double,1> sv_initial;
public:
  TcconSounding1Fixture();
  virtual ~TcconSounding1Fixture()
  {
    config_state_vector->update_state(sv_initial);
  }
  boost::shared_ptr<StateVector> config_state_vector;
  boost::shared_ptr<InitialGuess> config_initial_guess;
  boost::shared_ptr<Absorber> config_absorber;
  LuabindObject lua_config;
};

boost::shared_ptr<LuaState> TcconSounding1Fixture::ls;

TcconSounding1Fixture::TcconSounding1Fixture()
{ 
  if(!ls) {
    std::string ecmwf = test_data_dir() + "../test/tccon_small_set/acos_EcmB2900_tccon_5_good_qual.h5";
    std::string spectrum = test_data_dir() + "../test/tccon_small_set/acos_L1bB2900_tccon_5_good_qual.h5";
    std::string config = test_data_dir() + "../input/gosat/config/config.lua";
    setenv("met_file", ecmwf.c_str(), 1);
    setenv("spectrum_file", spectrum.c_str() , 1);
    setenv("sounding_id", "20100223034944", 1);
    ls = LuaState::load_file(config);
  }
  lua_config = ls->globals()["config"];
  config_state_vector = lua_config["state_vector"].value_ptr<StateVector>();
  config_initial_guess = 
    lua_config["initial_guess"].value_ptr<InitialGuess>();
  config_absorber = lua_config["absorber"].value_ptr<Absorber>();
  sv_initial.reference(config_initial_guess->initial_guess());
}

BOOST_FIXTURE_TEST_SUITE(absorber_absco_tccon1, TcconSounding1Fixture)
BOOST_AUTO_TEST_CASE(check_integration)
{
  AbsorberAbsco& a = dynamic_cast<AbsorberAbsco&>(*config_absorber);

  // From more extensive tests in
  // verification_test/gas_optical_depth_test/gas_optical_depth_test.py
  // we determined that the maximum difference with the integration is
  // at this value, for the top layer.
  double wnmax = 13056.32;
  int spec_index = 0;
  int species_index = 1;
  for(int layer = 0; layer < a.number_layer(); ++layer) {
    double od = 
      a.optical_depth_each_layer_nder(wnmax, spec_index)(layer,species_index);
    double od_direct =
      a.optical_depth_each_layer_direct_integrate(wnmax, spec_index, 
						  species_index, layer);
    // Agree to 0.2%
    BOOST_CHECK_CLOSE(od, od_direct, 0.25);
  }
}

BOOST_AUTO_TEST_CASE(timing)
{
  is_timing_test();
  AbsorberAbsco& a = dynamic_cast<AbsorberAbsco&>(*config_absorber);
  boost::timer tm;
  std::cerr << "Starting gas optical depth of all data\n";
  int j = 0;
  for(double i = 12929.94; i < 13210.15; i += 0.01, j++) {
    if(j % 1000 == 0)
      std::cerr << "Reading " << j << "\n"
		<< "Total time: " << tm.elapsed() << "\n";
    ArrayAd<double, 2> od = a.optical_depth_each_layer(i, 0);
  }
  std::cerr << "Done\n"
	    << "Total time: " << tm.elapsed() << "\n";
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(absorber_absco_conf, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(optical_depth_each_layer)
{
  secondIndex i2;
  IfstreamCs expected_data(test_data_dir() + 
		   "expected/absorber_absco/optical_depth_each_layer");
  Absorber& a = *config_absorber;
  Array<double, 1> od_expect;
  expected_data >> od_expect;
  Array<double, 1> od_calc_1( sum(a.optical_depth_each_layer(12929.94,0).value(),i2) );
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_1, od_expect, 1e-10);
  expected_data >> od_expect;
  Array<double, 1> od_calc_2( sum(a.optical_depth_each_layer(12930.30,0).value(),i2) );
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_2, od_expect, 1e-10);
  if (false) {
    std::cerr << setprecision(20) << std::scientific
	      << "# od_expect 12929.94" << std::endl
	      << od_calc_1 << std::endl
	      << "# od_expect 12930.30" << std::endl
	      << od_calc_2 << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(optical_depth_each_layer_direct_integrate)
{
  secondIndex i2;
  IfstreamCs expected_data(test_data_dir() + 
		   "expected/absorber_absco/optical_depth_each_layer");
  AbsorberAbsco& a = dynamic_cast<AbsorberAbsco&>(*config_absorber);
  Array<double, 1> od_expect;
  expected_data >> od_expect;
  Array<double, 1> od_calc_1( sum(a.optical_depth_each_layer_direct_integrate(12929.94,0),i2) );
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_1, od_expect, 1e-8);
  expected_data >> od_expect;
  Array<double, 1> od_calc_2( sum(a.optical_depth_each_layer_direct_integrate(12930.30,0),i2) );
  BOOST_CHECK_MATRIX_CLOSE_TOL(od_calc_2, od_expect, 1e-8);
}

BOOST_AUTO_TEST_CASE(press_wf)
{
  IfstreamCs expected_data(test_data_dir() + 
			   "expected/absorber_absco/pressure_wf");
  AbsorberAbsco& a = dynamic_cast<AbsorberAbsco&>(*config_absorber);
  Array<double, 1> press_wf_lay_expect, press_wf_lev_expect;
  expected_data >> press_wf_lay_expect >> press_wf_lev_expect;
  BOOST_CHECK_MATRIX_CLOSE(a.pressure_weighting_function_layer().value(), 
			   press_wf_lay_expect);
  BOOST_CHECK_MATRIX_CLOSE(a.pressure_weighting_function_grid().value(), 
			   press_wf_lev_expect);
  BOOST_CHECK_CLOSE(a.xgas("CO2").value(), 0.00037256930599168425, 1e-4);
  if (false) {
    std::cerr << setprecision(20) << std::scientific
	      << "# Expected pressure weighting function for layers." << std::endl
	      << a.pressure_weighting_function_layer().value() << std::endl
	      << "# Expected pressure weighting function for levels" << std::endl
	      << a.pressure_weighting_function_grid().value() << std::endl;
  }    
}

BOOST_AUTO_TEST_CASE(gas_column_thickness)
{
  // Make sure that wet - dry column == h2o column as it should
  AbsorberAbsco& a = dynamic_cast<AbsorberAbsco&>(*config_absorber);
  Array<double, 1> wet_col(a.wet_air_column_thickness_layer().value.value());
  Array<double, 1> dry_col(a.dry_air_column_thickness_layer().value.value());
  Array<double, 1> h2o_col(a.gas_column_thickness_layer("H2O").value.value());
  // The two comparison values are on the order of 1e23 so must have a higher tol here
  BOOST_CHECK_MATRIX_CLOSE_TOL(wet_col - dry_col, h2o_col, 1e13);

  // Check that dry / column(O2) = vmr
  Array<double, 1> o2_col(a.gas_column_thickness_layer("O2").value.value());
  Array<double, 1> o2_vmr(o2_col.rows());
  o2_vmr = a.average_vmr("O2").value();
  BOOST_CHECK_MATRIX_CLOSE(o2_col / dry_col, o2_vmr);
}

BOOST_AUTO_TEST_CASE(jacobian)
{
  AbsorberAbsco& a = dynamic_cast<AbsorberAbsco&>(*config_absorber);
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  // Pick an band where we have some CO2, so varying the VMR of CO2 affects 
  // the results.
  int spec_index = 2;
  double wn = 4820.0;
  ArrayAd<double, 2> od = a.optical_depth_each_layer(wn, spec_index);
  Array<double, 2> od0(od.shape());
  od0 = od.value();
  Array<double, 3> jac = od.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 2> jacfd(od0.shape());
    jacfd = (a.optical_depth_each_layer_nder(wn,spec_index) - od0) 
      / epsilon(i);
    if(false) {			// Can turn this off to dump values,
				// if needed for debugging
      double diff = max(abs(jac(Range::all(), Range::all(), i) - jacfd));
      if(diff > 0)
	std::cerr << i << ": " << jac(Range::all(), Range::all(), i) << "\n"
		  << jacfd << "\n"
		  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), Range::all(), i), jacfd, 
				 1e-10);
  }
}

BOOST_AUTO_TEST_SUITE_END()

