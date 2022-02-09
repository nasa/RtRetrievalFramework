#include "lidort_driver.h"
#include "unit_test_support.h"
#include "lidort_fixture.h"
#include "old_constant.h"

#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;

bool check_brdf_inputs(boost::shared_ptr<LidortRtDriver>& lidort_driver) {
  Lidort_Brdf_Sup_Accessories brdf_check = Lidort_Brdf_Sup_Accessories(
      lidort_driver->brdf_interface()->brdf_sup_in_ptr(),
      lidort_driver->lidort_interface()->lidort_fixin_ptr(),
      lidort_driver->lidort_interface()->lidort_modin_ptr());
  brdf_check.brdf_input_check();

  Lidort_Exception_Handling& brdf_check_status = brdf_check.lidort_brdfcheck_status();
  Lidort_Pars lid_pars = Lidort_Pars::instance();

  if (brdf_check_status.ts_status_inputcheck() != lid_pars.lidort_success)
    std::cerr << brdf_check_status << std::endl;

  return brdf_check_status.ts_status_inputcheck() == lid_pars.lidort_success;
}

BOOST_FIXTURE_TEST_SUITE(lidort_driver_lambertian_simple, LidortDriverLambertianFixture)

BOOST_AUTO_TEST_CASE(simple)
{
  int nlayer = 1;
  Array<double, 1> heights(nlayer+1);
  Array<double, 1> surface_params(1); 
  Array<double, 1> od(nlayer);
  Array<double, 1> ssa(nlayer);
  Array<double, 2> pf(lidort_driver->number_moment(),nlayer);
  double taug, taur;
  Range all = Range::all();

  double refl_calc;
  double refl_expected;

  // Turn off delta-m scaling
  lidort_driver->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  lidort_driver->set_plane_parallel();

  // Simple height grid evenly spaced
  heights(0) = 100;
  for(int hidx = 1; hidx < nlayer+1; hidx++) {
    heights(hidx) = heights(hidx-1) - heights(0)/nlayer;
  }

  // No aerosols, and depolization factor = 0 
  // so simplified phase function moments:
  pf = 0.0;
  pf(0, all) = 1.0;
  pf(2, all) = 0.5;

  // Use simple lambertian throughout
  int surface_type = LAMBERTIAN;
  surface_params = 0.0;

  ////////////////
  // Surface only
  surface_params(0) = 1.0;

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa = taur / od;

  refl_calc = lidort_driver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params,
                                                   od, ssa, pf);
  
  // Surface only = 1/pi
  // This checks % differerence, so tol is % diff
  refl_expected = 1/OldConstant::pi * surface_params(0);
  BOOST_CHECK_CLOSE(refl_expected, refl_calc, 1e-3);

  ////////////////
  // Gas + Surface
  surface_params(0) = 1.0;

  taur = 1.0e-6/nlayer;
  taug = 1.0/nlayer;

  od = taur + taug;
  ssa = taur / od;

  refl_calc = lidort_driver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params,
                                                   od, ssa, pf);

  refl_expected = 1/OldConstant::pi * exp(-1/cos(sza(0))) * exp(-1/cos(zen(0)));
  BOOST_CHECK_CLOSE(refl_expected, refl_calc, 1e-3);

  ////////////////
  // Rayleigh only
  surface_params(0) = 0.0;

  taur = 2.0e-2/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa = taur / od;

  refl_calc = lidort_driver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params,
                                                   od, ssa, pf);


  // Expected value from VLIDORT, this is
  // not going to agree to infinite precision
  // Note the tolerance here is in %
  refl_expected = 2.387246757232095E-003;
  BOOST_CHECK_CLOSE(refl_expected, refl_calc, 6e-3);
}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE(lidort_driver_coxmunk_simple, LidortDriverCoxmunkFixture)

BOOST_AUTO_TEST_CASE(simple)
{
  int nlayer = 1;
  Array<double, 1> heights(nlayer+1);
  Array<double, 1> surface_params(5); 
  Array<double, 1> od(nlayer);
  Array<double, 1> ssa(nlayer);
  Array<double, 2> pf(lidort_driver->number_moment(), nlayer);
  double taug, taur;
  Range all = Range::all();

  double refl_calc;
  double refl_expected;

  // Turn off delta-m scaling
  lidort_driver->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  lidort_driver->set_plane_parallel();

  // Simple height grid evenly spaced
  heights(0) = 100;
  for(int hidx = 1; hidx < nlayer+1; hidx++) {
    heights(hidx) = heights(hidx-1) - heights(0)/nlayer;
  }

  // No aerosols, and depolization factor = 0 
  // so simplified phase function moments:
  pf = 0.0;
  pf(0, all) = 1.0;
  pf(2, all) = 0.5;

  // Use simple lambertian throughout
  int surface_type = COXMUNK;
  surface_params = 0.0;

  ////////////////
  // Surface only
  surface_params(0) = 1.0;
  surface_params(1) = 1.0e-6;
  surface_params(2) = 1.334;
  
  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa = taur / od;

  refl_calc = lidort_driver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params,
                                                   od, ssa, pf);

  BOOST_CHECK_EQUAL(check_brdf_inputs(lidort_driver), true);
  
  // Value for VLIDORT
  refl_expected = 0.54319850628416033;
  BOOST_CHECK_CLOSE(refl_expected, refl_calc, 1e-3);

}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE(lidort_driver_coxmunk_plus_lamb_simple, LidortDriverCoxmunkFixture)

BOOST_AUTO_TEST_CASE(simple)
{
  int nlayer = 1;
  Array<double, 1> heights(nlayer+1);
  Array<double, 1> surface_params(5); 
  ArrayAd<double, 1> od(nlayer, 1);
  ArrayAd<double, 1> ssa(nlayer, 1);
  ArrayAd<double, 2> pf(lidort_driver->number_moment(),nlayer, 1);
  double taug, taur;
  Range all = Range::all();

  double refl_calc;
  double refl_expt;

  // Simple height grid evenly spaced
  heights(0) = 100;
  for(int hidx = 1; hidx < nlayer+1; hidx++) {
    heights(hidx) = heights(hidx-1) - heights(0)/nlayer;
  }

  // No aerosols, and depolization factor = 0 
  // so simplified phase function moments:
  pf = 0.0;
  pf(0, all) = 1.0;
  pf(2, all) = 0.5;

  // Use simple lambertian throughout
  int surface_type = COXMUNK;
  surface_params = 0.0;

  ////////////////
  // Surface only
  surface_params(0) = 1.0;
  surface_params(1) = 1.0e-6;
  surface_params(2) = 1.334;
  surface_params(3) = 0.5;

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa = 0;
  ssa.value() = taur / od.value();
  
  //////
  // Plane-parallel
  lidort_driver->set_plane_parallel();

  blitz::Array<double, 2> jac_atm;
  blitz::Array<double, 1> jac_surf;
  ArrayAd<double, 1> lidort_surface(surface_params.shape(), 1);

  lidort_surface.value() = surface_params;
  lidort_surface.jacobian() = 1.0;

  lidort_driver->reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                    surface_type, lidort_surface,
                                                    od, ssa, pf, refl_calc, jac_atm, jac_surf);

  BOOST_CHECK_EQUAL(check_brdf_inputs(lidort_driver), true);

  // Adjust analytic jacobians have same meaning as fd jacobians
  jac_surf(1) *= lidort_surface.jacobian()(1,0);
  jac_surf(2) *= lidort_surface.jacobian()(2,0);

  // Value for VLIDORT
  refl_expt = 0.70235315460259928;
  BOOST_CHECK_CLOSE(refl_expt, refl_calc, 1e-3);

  // Check surface jacobians against FD

  blitz::Array<double, 1> pert_values(surface_params.extent(firstDim)-1);
  pert_values = 0;
  pert_values = 1e-8, 1e-8, 1e-8, 1e-6;

  blitz::Array<double, 1> jac_surf_fd( jac_surf.extent() );
  double refl_fd;

  jac_surf_fd = 0.0;
  refl_fd = 0.0;

  // First check PP mode against value just computed
  for(int p_idx = 0; p_idx < pert_values.extent(firstDim); p_idx++) {
    blitz::Array<double,1> surface_params_pert( surface_params.extent() );
    surface_params_pert = surface_params;
    surface_params_pert(p_idx) += pert_values(p_idx);

    refl_fd = lidort_driver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params_pert,
                                                   od.value(), ssa.value(), pf.value());

    jac_surf_fd(p_idx) = (refl_fd - refl_calc) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf, jac_surf_fd, 1e-7);

  // Pseudo-spherical mode FD test
  lidort_driver->set_pseudo_spherical();
 
  lidort_surface.value() = surface_params;
  lidort_surface.jacobian() = 1.0;

  lidort_driver->reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                    surface_type, lidort_surface,
                                                    od, ssa, pf, refl_calc, jac_atm, jac_surf);

  BOOST_CHECK_EQUAL(check_brdf_inputs(lidort_driver), true);

  // Adjust analytic jacobians have same meaning as fd jacobians
  jac_surf(1) *= lidort_surface.jacobian()(1,0);
  jac_surf(2) *= lidort_surface.jacobian()(2,0);

  for(int p_idx = 0; p_idx < pert_values.extent(firstDim); p_idx++) {
    blitz::Array<double,1> surface_params_pert( surface_params.extent() );
    surface_params_pert = surface_params;
    surface_params_pert(p_idx) += pert_values(p_idx);

    refl_fd = lidort_driver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params_pert,
                                                   od.value(), ssa.value(), pf.value());

    jac_surf_fd(p_idx) = (refl_fd - refl_calc) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf, jac_surf_fd, 1e-7);

  // Line-of-site mode FD test
  sza = sza + 1e-3; // or else risk divide by zero
  azm = azm + 1e-3; // or else risk divide by zero
  zen = zen + 1e-3; // or else risk divide by zero

  lidort_driver->set_line_of_sight();

  lidort_surface.value() = surface_params;
  lidort_surface.jacobian() = 1.0;
  lidort_driver->reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                 surface_type, lidort_surface,
                                                 od, ssa, pf, refl_calc, jac_atm, jac_surf);

  BOOST_CHECK_EQUAL(check_brdf_inputs(lidort_driver), true);

  // Adjust analytic jacobians have same meaning as fd jacobians
  jac_surf(1) *= lidort_surface.jacobian()(1,0);
  jac_surf(2) *= lidort_surface.jacobian()(2,0);

  for(int p_idx = 0; p_idx < pert_values.extent(firstDim); p_idx++) {
    blitz::Array<double,1> surface_params_pert( surface_params.extent() );
    surface_params_pert = surface_params;
    surface_params_pert(p_idx) += pert_values(p_idx);

    refl_fd = lidort_driver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params_pert,
                                                   od.value(), ssa.value(), pf.value());

    jac_surf_fd(p_idx) = (refl_fd - refl_calc) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf, jac_surf_fd, 1e-7);

}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(lidort_driver_brdf_veg, LidortDriverCommonFixture)

BOOST_AUTO_TEST_CASE(simple)
{

  int nstreams = 4;
  int nmoms = 2*nstreams;
  bool do_multiple_scattering_only = false;

  // Use simple BREON throughout
  int surface_type = BREONVEG;
  ArrayAd<double, 1> surface_params(5, 1); 
  surface_params = 0.0;

  boost::shared_ptr<LidortRtDriver> lidort_driver(new LidortRtDriver(nstreams, nmoms, do_multiple_scattering_only, surface_type, zen, pure_nadir));

  int nlayer = 1;
  Array<double, 1> heights(nlayer+1);
  ArrayAd<double, 1> od(nlayer,1);
  ArrayAd<double, 1> ssa(nlayer,1);
  ArrayAd<double, 2> pf(lidort_driver->number_moment(),nlayer,1);
  double taug, taur;
  Range all = Range::all();

  double refl_calc;
  blitz::Array<double, 2> jac_atm;
  blitz::Array<double, 1> jac_surf;

  // Turn off delta-m scaling
  lidort_driver->lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  lidort_driver->set_plane_parallel();

  // Simple height grid evenly spaced
  heights(0) = 100;
  for(int hidx = 1; hidx < nlayer+1; hidx++) {
    heights(hidx) = heights(hidx-1) - heights(0)/nlayer;
  }

  // No aerosols, and depolization factor = 0 
  // so simplified phase function moments:
  pf = 0.0;
  pf(0, all) = 1.0;
  pf(2, all) = 0.5;

  ////////////////
  // Surface only
  surface_params(0) = 1.0; // rahman kernel factor
  surface_params(1) = 0.1; // hotspot parameter
  surface_params(2) = 0.3; // asymmetry
  surface_params(3) = 1.5; // anisotropy_parameter
  surface_params(4) = 1.0; // breon kernel factor

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  od = taur + taug;
  ssa = 0;
  ssa.value() = taur / od.value();

  // Set sza to compare with that used in the l_rad test we compare against
  sza = 0.1;

  lidort_driver->reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                    surface_type, surface_params,
                                                    od, ssa, pf, refl_calc, jac_atm, jac_surf);

  // Compare against an offline calculated value, or could compare against value from l_rad
  double refl_expected = 0.035435854422713485;
  BOOST_CHECK_CLOSE(refl_expected, refl_calc, 1e-3);

  // Check surface jacobians against FD

  blitz::Array<double, 1> pert_values(surface_params.rows());
  pert_values = 1e-8, 1e-8, 1e-6, 1e-6, 1e-6;

  blitz::Array<double, 1> jac_surf_fd( jac_surf.extent() );
  double refl_fd;

  jac_surf_fd = 0.0;
  refl_fd = 0.0;

  // First check PP mode against value just computed
  for(int p_idx = 0; p_idx < pert_values.extent(firstDim); p_idx++) {
    blitz::Array<double,1> surface_params_pert( surface_params.rows() );
    surface_params_pert = surface_params.value();
    surface_params_pert(p_idx) += pert_values(p_idx);

    refl_fd = lidort_driver->reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params_pert,
                                                   od.value(), ssa.value(), pf.value());

    jac_surf_fd(p_idx) = (refl_fd - refl_calc) / pert_values(p_idx);
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf, jac_surf_fd, 2e-7);

}

BOOST_AUTO_TEST_SUITE_END()
