#include "unit_test_support.h"
#include "atmosphere_fixture.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(atmosphere_oco, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  IfstreamCs expected_data(test_data_dir() + 
			   "expected/atmosphere_oco/rt_parameters_each_layer");
  // Expected values were gotten by running the old Fortran code and
  // extracting out the answer from that.
  Array<double, 1> od_expect, ssa_expect;
  // Pick one level to check, because the matrix is big. No
  // significance in the level we selected here, it was just one in
  // the middle.
  Array<double, 2> scat_momsub_expect;
  expected_data >> scat_momsub_expect;
  expected_data >> od_expect;
  expected_data >> ssa_expect;
  ArrayAd<double, 2> iv(atm->intermediate_variable(12929.94, 0));
  BOOST_CHECK_MATRIX_CLOSE_TOL(atm->optical_depth_wrt_iv(12929.94, 0).value(), 
			       od_expect, 1e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL(atm->optical_depth_wrt_iv(12929.94, 0,iv).value(), 
			       od_expect, 1e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (atm->single_scattering_albedo_wrt_iv(12929.94, 0).value(), 
     ssa_expect, 1e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL
    (atm->single_scattering_albedo_wrt_iv(12929.94, 0,iv).value(), ssa_expect, 
     1e-6);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12929.94, 0, "O2").value(), 
		    0.00043001078262954612, 1e-4);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12929.94, 0, "H2O").value(),
		    0, 1e-4);
  Array<double, 2> scat_momsub
    (atm->scattering_moment_wrt_iv(12929.94, 0, 4, 1).value()
     (Range::all(), 9, Range::all()));
  BOOST_CHECK_MATRIX_CLOSE_TOL(scat_momsub, 
			       scat_momsub_expect, 1e-4);
  scat_momsub = (atm->scattering_moment_wrt_iv(12929.94, 0, iv, 4, 1).value()
		 (Range::all(), 9, Range::all()));
  BOOST_CHECK_MATRIX_CLOSE_TOL(scat_momsub, 
			       scat_momsub_expect, 1e-4);
  expected_data >> scat_momsub_expect;
  expected_data >> od_expect;
  expected_data >> ssa_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(atm->optical_depth_wrt_iv(12930.30, 0).value(), 
			       od_expect, 1e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL(atm->single_scattering_albedo_wrt_iv(12930.30, 0).value(),
			       ssa_expect, 1e-6);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12930.30, 0, "O2").value(), 
		    0.00043415109976773604, 1e-4);
  BOOST_CHECK_CLOSE(atm->column_optical_depth(12930.30, 0, "H2O").value(),
		    0, 1e-4);
  Array<double, 2> scat_momsub2(atm->scattering_moment_wrt_iv(12930.30, 0, 4, 1).value()
				(Range::all(), 9, Range::all()));
  BOOST_CHECK_MATRIX_CLOSE_TOL(scat_momsub2, 
			       scat_momsub_expect, 1e-4);
  BOOST_CHECK_EQUAL(count(statev->used_flag()), 109 - 5);
}

BOOST_AUTO_TEST_CASE(rayleigh_atmosphere)
{
  // We check that leaving out the aerosol gives the same results as
  // simply having the aerosol extinction set to 0.
  boost::shared_ptr<AtmosphereOco> atm_zeroext = atm->clone();
  // Set state vector so that the extinction coefficient of
  // atm_zeroext is 0.
  StateVector sv;
  sv.add_observer(*atm_zeroext->aerosol_ptr());
  Array<double, 1> x(sv.observer_claimed_size());
  x = -1000;
  sv.update_state(x);

  boost::shared_ptr<Pressure> pressure_clone = atm->pressure_ptr()->clone();
  boost::shared_ptr<Temperature> temperature_clone = 
    atm->temperature_ptr()->clone(pressure_clone);
  boost::shared_ptr<Ground> ground_clone = atm->ground()->clone();
  std::vector<boost::shared_ptr<Altitude> > alt_clone;
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, atm->altitude_ptr())
    alt_clone.push_back(a->clone(pressure_clone, temperature_clone));
  boost::shared_ptr<Absorber> absorber_clone =
    atm->absorber_ptr()->clone(pressure_clone, temperature_clone, alt_clone);
  boost::shared_ptr<RelativeHumidity> rh_clone =
    atm->relative_humidity_ptr()->clone(absorber_clone, temperature_clone, pressure_clone);
  boost::shared_ptr<AerosolOptical> aerosol_null;
  boost::shared_ptr<AtmosphereOco> 
    atm_rayleigh(new AtmosphereOco(absorber_clone,
				   pressure_clone,
				   temperature_clone,
				   aerosol_null,
				   rh_clone,
				   ground_clone,
				   alt_clone,
				   atm->constant_ptr()));
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->ground()->surface_parameter(12929.94, 0).value(),
			   atm_zeroext->ground()->surface_parameter(12929.94, 0).value());
  BOOST_CHECK_EQUAL(atm_rayleigh->number_spectrometer(),
		    atm_zeroext->number_spectrometer());
  BOOST_CHECK_EQUAL(atm_rayleigh->number_layer(),
		    atm_zeroext->number_layer());
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->optical_depth_wrt_iv(12929.94, 0).value(),
			   atm_zeroext->optical_depth_wrt_iv(12929.94, 0).value());
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->single_scattering_albedo_wrt_iv(12929.94, 0).value(), 
			   atm_zeroext->single_scattering_albedo_wrt_iv(12929.94, 0).value());
  BOOST_CHECK_CLOSE
    (atm_rayleigh->column_optical_depth(12929.94, 0, "O2").value(), 
     atm_zeroext->column_optical_depth(12929.94, 0, "O2").value(), 
     1e-6);
  BOOST_CHECK_CLOSE
    (atm_rayleigh->column_optical_depth(12929.94, 0, "H2O").value(),
     atm_zeroext->column_optical_depth(12929.94, 0, "H2O").value(),
     1e-6);
  Range r1(0,atm_rayleigh->scattering_moment_wrt_iv(12738.853381475927, 0).rows() - 1);
  Range r2(0,
	   atm_rayleigh->scattering_moment_wrt_iv(12738.853381475927, 0).depth() - 1);
  Range ra = Range::all();
  BOOST_CHECK_MATRIX_CLOSE(atm_rayleigh->scattering_moment_wrt_iv(12929.94, 0).value(),
			   atm_zeroext->scattering_moment_wrt_iv(12929.94, 0).value()(r1, ra, r2));
}

BOOST_AUTO_TEST_CASE(uplooking_atmosphere)
{
  // We check that leaving out the ground gives the same results 

  boost::shared_ptr<Pressure> pressure_clone = atm->pressure_ptr()->clone();
  boost::shared_ptr<Temperature> temperature_clone = 
    atm->temperature_ptr()->clone(pressure_clone);
  boost::shared_ptr<Aerosol> aerosol_clone =
    atm->aerosol_ptr()->clone();
  std::vector<boost::shared_ptr<Altitude> > alt_clone;
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, atm->altitude_ptr())
    alt_clone.push_back(a->clone(pressure_clone, temperature_clone));
  boost::shared_ptr<Absorber> absorber_clone =
    atm->absorber_ptr()->clone(pressure_clone, temperature_clone, alt_clone);
  boost::shared_ptr<RelativeHumidity> rh_clone =
    atm->relative_humidity_ptr()->clone(absorber_clone, temperature_clone, pressure_clone);
  boost::shared_ptr<Ground> ground_null;
  boost::shared_ptr<AtmosphereOco> 
    atm_uplooking(new AtmosphereOco(absorber_clone,
				    pressure_clone,
				    temperature_clone,
				    aerosol_clone,
				    rh_clone,
				    ground_null,
				    alt_clone,
				    atm->constant_ptr()));

  BOOST_CHECK_EQUAL(atm_uplooking->number_spectrometer(),
		    atm->number_spectrometer());
  BOOST_CHECK_EQUAL(atm_uplooking->number_layer(),
		    atm->number_layer());
  BOOST_CHECK_MATRIX_CLOSE(atm_uplooking->optical_depth_wrt_iv(12929.94, 0).value(),
			   atm->optical_depth_wrt_iv(12929.94, 0).value());
  BOOST_CHECK_MATRIX_CLOSE(atm_uplooking->single_scattering_albedo_wrt_iv(12929.94, 0).value(),
			   atm->single_scattering_albedo_wrt_iv(12929.94, 0).value());
  BOOST_CHECK_CLOSE
    (atm_uplooking->column_optical_depth(12929.94, 0, "O2").value(), 
     atm->column_optical_depth(12929.94, 0, "O2").value(), 1e-6);
  BOOST_CHECK_CLOSE
    (atm_uplooking->column_optical_depth(12929.94, 0, "H2O").value(),
     atm->column_optical_depth(12929.94, 0, "H2O").value(),
     1e-6);
  BOOST_CHECK_MATRIX_CLOSE(atm_uplooking->scattering_moment_wrt_iv(12929.94, 0).value(),
			   atm->scattering_moment_wrt_iv(12929.94, 0).value());
  BOOST_CHECK_MATRIX_CLOSE(atm_uplooking->scattering_moment_wrt_iv(12929.94, 0).jacobian(),
			   atm->scattering_moment_wrt_iv(12929.94, 0).jacobian());
  Array<double, 2> l_opdel_uplooking = 
    atm_uplooking->optical_depth_wrt_iv(12929.94, 0).jacobian();
  Array<double, 2> l_opdel =
    atm->optical_depth_wrt_iv(12929.94, 0).jacobian();
  Array<double, 2> l_ssa_uplooking = 
    atm_uplooking->single_scattering_albedo_wrt_iv(12929.94, 0).jacobian();
  Array<double, 2> l_ssa =
    atm->single_scattering_albedo_wrt_iv(12929.94, 0).jacobian();
  BOOST_CHECK_MATRIX_CLOSE(l_opdel_uplooking, l_opdel);
  BOOST_CHECK_MATRIX_CLOSE(l_ssa_uplooking, l_ssa);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(atmosphere_oco_jac, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(optical_depth_jac)
{
  is_long_test();
  RtAtmosphere& atm = *config_atmosphere;
  IfstreamCs expected_data(test_data_dir() + 
			   "expected/atmosphere_oco/rt_parameters_each_layer");
  Array<double, 1> od_expect, ssa_expect;
  Array<double, 2> scat_momsub_expect;
  expected_data >> scat_momsub_expect;
  expected_data >> od_expect;
  expected_data >> ssa_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(
	   atm.optical_depth_wrt_state_vector(12929.94, 0).value(), 
	   od_expect, 1e-6);
  // Pick an band where we have some CO2, so varying the VMR of CO2 affects 
  // the results.
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  int spec_index = 2;
  double wn = 4820.0;
  ArrayAd<double, 1> od = 
    atm.optical_depth_wrt_state_vector(wn, spec_index);
  Array<double, 1> od0(od.shape());
  od0 = od.value();
  Array<double, 2> jac = od.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 1> jacfd(od0.shape());
    jacfd = (atm.optical_depth_wrt_iv(wn,spec_index).value() - od0) 
      / epsilon(i);
    if(false) {			// Can turn this off to dump values,
				// if needed for debugging
      double diff = max(abs(jac(Range::all(), i) - jacfd));
      if(diff > 0)
	std::cerr << i << ": " << jac(Range::all(), i) << "\n"
		  << jacfd << "\n"
		  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), i), jacfd, 
				 5e-9);
  }
}

BOOST_AUTO_TEST_CASE(ssa_jac)
{
  is_long_test();
  RtAtmosphere& atm = *config_atmosphere;
  IfstreamCs expected_data(test_data_dir() + 
			   "expected/atmosphere_oco/rt_parameters_each_layer");
  Array<double, 1> od_expect, ssa_expect;
  Array<double, 2> scat_momsub_expect;
  expected_data >> scat_momsub_expect;
  expected_data >> od_expect;
  expected_data >> ssa_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(
    atm.single_scattering_albedo_wrt_state_vector(12929.94, 0).value(), 
    ssa_expect, 1e-6);
  // Pick an band where we have some CO2, so varying the VMR of CO2 affects 
  // the results.
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  int spec_index = 2;
  double wn = 4820.0;
  ArrayAd<double, 1> ssa = 
    atm.single_scattering_albedo_wrt_state_vector(wn, spec_index);
  Array<double, 1> ssa0(ssa.shape());
  ssa0 = ssa.value();
  Array<double, 2> jac = ssa.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 1> jacfd(ssa0.shape());
    jacfd = (atm.single_scattering_albedo_wrt_iv(wn,spec_index).value() - ssa0) 
      / epsilon(i);
    if(false) {			// Can turn this off to dump values,
				// if needed for debugging
      double diff = max(abs(jac(Range::all(), i) - jacfd));
      if(diff > 0)
	std::cerr << i << ": " << jac(Range::all(), i) << "\n"
		  << jacfd << "\n"
		  << diff << "\n";
    }
    // There are a wide range in the size of the difference, because
    // the Jacobian values vary wildly in size from one row to the
    // next. We looked at all the values, and things look
    // reasonable. To have an automated test, we look at having either
    // the difference being very small, or the relative difference
    // small.
    Array<double, 1> diff(jac(Range::all(), i) - jacfd);
    BOOST_CHECK(max(abs(diff)) < 5e-8 ||
		max(abs(where(jacfd == 0, 0, diff / jacfd))) < 5e-4);
  }
}

BOOST_AUTO_TEST_CASE(scattering_moment_jac)
{
  is_long_test();
  RtAtmosphere& atm = *config_atmosphere;
  StateVector& sv = *config_state_vector;
  Array<double, 1> sv0(sv.state().copy());
  int spec_index = 2;
  double wn = 4820.0;
  ArrayAd<double, 3> sm = 
    atm.scattering_moment_wrt_state_vector(wn, spec_index);
  Array<double, 3> sm0(sm.shape());
  sm0 = sm.value();
  Array<double, 4> jac = sm.jacobian().copy();
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 3> jacfd(sm.shape());
    jacfd = (atm.scattering_moment_wrt_iv(wn,spec_index).value() - sm0) 
      / epsilon(i);
    Array<double, 3> diff(jac(Range::all(), Range::all(), Range::all(), i) - 
			  jacfd);
    if(false) {			// Can turn this off to dump values,
				// if needed for debugging
      if(max(abs(diff)) > 0) {
	std::cerr << i << ": " << max(abs(diff)) << " " 
		  << max(abs(where(abs(jacfd) < 1e-15, 0, diff / jacfd))) 
		  << "\n";
	std::cerr << jac(Range(0,9), Range::all(), 0, i) << "\n"
		  << jacfd(Range(0,9), Range::all(), 0) << "\n";
      }
    }

  BOOST_CHECK(all(abs(diff) < 2e-9 || 
		  abs(where(jacfd == 0, 0, diff / jacfd)) < 1e-4));
  }
}

BOOST_AUTO_TEST_CASE(optical_depth_timing)
{
  is_timing_test();
  RtAtmosphere& atm = *config_atmosphere;
  int i = 0;
  for(double wn = 12929.94; wn <= 13210.15; wn += 0.01) {
    ArrayAd<double, 1> od = atm.optical_depth_wrt_iv(wn, 0);
    if(++i % 1000 == 0)
      std::cerr << "Done with " << i << "\n"
		<< atm.timer_info() << "\n";
  }
}

BOOST_AUTO_TEST_SUITE_END()
