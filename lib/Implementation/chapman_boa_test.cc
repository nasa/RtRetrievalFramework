#include "chapman_boa.h"
#include "unit_test_support.h"
#include "array_ad.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(chapman_boa, GlobalFixture)

BOOST_AUTO_TEST_CASE(tester)
{

  int nlayers = 23;

  // Get the preprepared atmosphere
  Array<double, 2> raymoms(3,nlayers);
  Array<double, 1> height_grid(nlayers+1);
  Array<double, 1> molext(nlayers);
  Array<double, 1> molomg(nlayers);

  IfstreamCs atm_data(test_data_dir() + "expected/chapman_boa/iop_dump5_gas_csky.res_l");
  int int_dummy;
  double dbl_dummy;
  for (int mom_idx = 0; mom_idx < 3; mom_idx++) {
    atm_data >> int_dummy;
    for(int lay_idx = 0; lay_idx < nlayers; lay_idx++)
      atm_data >> raymoms(mom_idx, lay_idx);
  }

  height_grid(0) = 60.0;
  for(int lay_idx = 0; lay_idx < nlayers; lay_idx++)
    atm_data >> int_dummy 
	     >> height_grid(lay_idx+1) 
	     >> molext(lay_idx) 
	     >> molomg(lay_idx)
	     >> dbl_dummy >> dbl_dummy >> dbl_dummy >> dbl_dummy;

  Array<double, 1> temperature_grid(nlayers+1);
  temperature_grid(0) = 220.0e0;
  for(int n = 1; n <= nlayers; n++)
    temperature_grid(n) = 220.0e0 + 3.5e0*double(n-1);
  
  Array<double, 1> pressure_grid(nlayers+1);
  pressure_grid(nlayers) = 1013.25e0;
  double scaleheight = log(2.0e0) / 5.5e0;
  for(int n = nlayers; n >= 1; n--) {
    double diff = height_grid(n-1) - height_grid(n);
    pressure_grid(n-1) = pressure_grid(n) * exp( - scaleheight * diff );
  }

  
  int n_szangles = 10;
  Array<double, 1> boa_szangles(n_szangles);
  boa_szangles = 40.0e0, 82.0e0, 85.0e0, 87.0e0, 87.5e0, 88.0e0, 88.5e0, 89.0e0, 89.5e0, 89.9e0;

  double rearth = 6371.0e0;
  double rfindex_parameter = 0.000288;

  Array<AutoDerivative<double>, 1> height_grid_ad(height_grid.extent(firstDim));
  height_grid_ad = cast<AutoDerivative<double> >(height_grid);
  Array<AutoDerivative<double>, 1> pressure_grid_ad(pressure_grid.extent(firstDim));
  pressure_grid_ad = cast<AutoDerivative<double> >(pressure_grid);
  pressure_grid_ad *= 1e2; // convert mb to Pa
  Array<AutoDerivative<double>, 1> temperature_grid_ad(temperature_grid.extent(firstDim));
  temperature_grid_ad = cast<AutoDerivative<double> >(temperature_grid);
  Array<AutoDerivative<double>, 1> molext_ad(molext.extent(firstDim));
  molext_ad = cast<AutoDerivative<double> >(molext);

  ////////
  // Check refractive geometry with pseudo-spherical correction
  std::vector<boost::shared_ptr<AtmRefractiveIndex> > refr_index;
  for(int beam_idx = 0; beam_idx < n_szangles; beam_idx++)
    refr_index.push_back(boost::shared_ptr<AtmRefractiveIndex>( new SimpleRefractiveIndex(rfindex_parameter, pressure_grid_ad, temperature_grid_ad)));
  ChapmanBOA boa_refractive = ChapmanBOA(rearth, boa_szangles, height_grid_ad, refr_index);

  Array<double, 1> boa_sza_expected;
  Array<double, 1> toa_sza_expected;
  Array<double, 1> entry_sza_expected;
  Array<double, 1> height_expected;
  Array<double, 1> extinction_expected;
  Array<double, 2> chapman_factors_expected;
  Array<double, 1> transmittances_expected;

  IfstreamCs expected_data_ps_refr(test_data_dir() + "expected/chapman_boa/offline_ps_refr.out");
  expected_data_ps_refr >> boa_sza_expected
			>> toa_sza_expected
			>> entry_sza_expected
			>> height_expected
			>> extinction_expected
			>> chapman_factors_expected
			>> transmittances_expected;

  // Check inputs
  BOOST_CHECK_MATRIX_CLOSE(boa_sza_expected, boa_szangles);
  BOOST_CHECK_MATRIX_CLOSE_TOL(height_expected(Range(0,nlayers-1)), height_grid(Range(1,nlayers)), 1e-5);

  // Check outputs
  BOOST_CHECK_MATRIX_CLOSE_TOL(chapman_factors_expected, value(boa_refractive.chapman_factors()), 1e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL(toa_sza_expected, value(boa_refractive.toa_nadir_szangles()), 1e-7);
  BOOST_CHECK_MATRIX_CLOSE_TOL(entry_sza_expected, value(boa_refractive.toa_entry_szangles()), 1e-8);
  BOOST_CHECK_MATRIX_CLOSE_TOL(extinction_expected, molext, 1e-8);
  BOOST_CHECK_MATRIX_CLOSE_TOL(transmittances_expected, value(boa_refractive.transmittance(molext_ad)), 1e-8);

  ////////
  // Check straight line geometry plane-parallel
  bool do_plane_parallel = true;
  ChapmanBOA boa_straight_pp = ChapmanBOA(rearth, do_plane_parallel, boa_szangles, height_grid_ad);

  IfstreamCs expected_data_pp_straight(test_data_dir() + "expected/chapman_boa/offline_pp_straight.out");
  expected_data_pp_straight >> boa_sza_expected
			    >> toa_sza_expected
			    >> entry_sza_expected
			    >> height_expected
			    >> extinction_expected
			    >> chapman_factors_expected
			    >> transmittances_expected;
  
  // Check inputs
  BOOST_CHECK_MATRIX_CLOSE(boa_sza_expected, boa_szangles);
  BOOST_CHECK_MATRIX_CLOSE_TOL(height_expected(Range(0,nlayers-1)), height_grid(Range(1,nlayers)), 1e-5);

  // Check outputs
  BOOST_CHECK_MATRIX_CLOSE_TOL(chapman_factors_expected, value(boa_straight_pp.chapman_factors()), 1e-5);
  BOOST_CHECK_MATRIX_CLOSE_TOL(toa_sza_expected, value(boa_straight_pp.toa_nadir_szangles()), 1e-5);
  BOOST_CHECK_MATRIX_CLOSE_TOL(entry_sza_expected, value(boa_straight_pp.toa_entry_szangles()), 1e-5);
  BOOST_CHECK_MATRIX_CLOSE_TOL(extinction_expected, molext, 1e-5);
  BOOST_CHECK_MATRIX_CLOSE_TOL(transmittances_expected, value(boa_straight_pp.transmittance(molext_ad)), 1e-4);

  ////////
  // Check straight line geometry pseudo-spherical
  do_plane_parallel = false;
  ChapmanBOA boa_straight_ps = ChapmanBOA(rearth, do_plane_parallel, boa_szangles, height_grid_ad);

  IfstreamCs expected_data_ps_straight(test_data_dir() + "expected/chapman_boa/offline_ps_straight.out");
  expected_data_ps_straight >> boa_sza_expected
			    >> toa_sza_expected
			    >> entry_sza_expected
			    >> height_expected
			    >> extinction_expected
			    >> chapman_factors_expected
			    >> transmittances_expected;
  
  // Check inputs
  BOOST_CHECK_MATRIX_CLOSE(boa_sza_expected, boa_szangles);
  BOOST_CHECK_MATRIX_CLOSE_TOL(height_expected(Range(0,nlayers-1)), height_grid(Range(1,nlayers)), 1e-5);

  // Check outputs
  BOOST_CHECK_MATRIX_CLOSE_TOL(chapman_factors_expected, value(boa_straight_ps.chapman_factors()), 1e-5);
  BOOST_CHECK_MATRIX_CLOSE_TOL(toa_sza_expected, value(boa_straight_ps.toa_nadir_szangles()), 1e-5);
  BOOST_CHECK_MATRIX_CLOSE_TOL(entry_sza_expected, value(boa_straight_ps.toa_entry_szangles()), 1e-5);
  BOOST_CHECK_MATRIX_CLOSE_TOL(extinction_expected, molext, 1e-5);
  BOOST_CHECK_MATRIX_CLOSE_TOL(transmittances_expected, value(boa_straight_ps.transmittance(molext_ad)), 1e-4);

}

BOOST_AUTO_TEST_SUITE_END()
