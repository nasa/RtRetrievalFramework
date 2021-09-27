#include "twostream_driver.h"
#include "unit_test_support.h"
#include "ground.h"
#include "old_constant.h"
#include "lidort_driver.h"
#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(twostream_driver, GlobalFixture)

void test_twostream(int surface_type, ArrayAd<double, 1>& surface_params, ArrayAd<double, 1>& taug, ArrayAd<double, 1>& taur, Array<double, 1>& pert_atm, Array<double, 1>& pert_surf, bool debug_output)
{
  blitz::Array<double, 1> sza, zen, azm;
  sza.resize(1); zen.resize(1); azm.resize(1);
  sza = 0.001;
  zen = 0.0;
  azm = 0.0;
  bool pure_nadir = false;

  int nlayer = taug.rows();
  int nparam = taug.number_variable();

  // Set up convenience ranges
  Range all = Range::all();
  Range rjac(0,nparam-1);
  Range rlay(0,nlayer-1);
  Range rsurf(0,surface_params.number_variable()-1);
  double refl_ts;
  double refl_lid;

  blitz::Array<double, 2> jac_atm_ts;
  blitz::Array<double, 1> jac_surf_ts;

  blitz::Array<double, 2> jac_atm_lid;
  blitz::Array<double, 1> jac_surf_lid;

  TwostreamRtDriver twostream_driver = TwostreamRtDriver(nlayer, surface_type, false);

  // Turn off delta-m scaling
  twostream_driver.twostream_interface()->do_d2s_scaling(false);

  // Plane-parallel
  twostream_driver.twostream_interface()->do_plane_parallel(true);

  // Use LIDORT for comparison
  int lid_nstreams = 1;
  int lid_nmoms = 3;
  bool do_multiple_scattering_only = true;
  LidortRtDriver lidort_driver = LidortRtDriver(lid_nstreams, lid_nmoms, 
                                                do_multiple_scattering_only,
                                                surface_type, zen, pure_nadir);  

  // Turn off delta-m scaling
  lidort_driver.lidort_interface()->lidort_modin().mbool().ts_do_deltam_scaling(false);

  // Plane-parallel
  lidort_driver.set_plane_parallel();

  // Simple height grid evenly spaced
  Array<double, 1> heights(nlayer+1);
  heights(0) = 100;
  for(int hidx = 1; hidx < nlayer+1; hidx++) {
    heights(hidx) = heights(hidx-1) - heights(0)/nlayer;
  }

  // Set up atmosphere jacobians to be one for taug and the other for taur
  taug.jacobian()(all, 0) = 1;
  taur.jacobian()(all, 1) = 1;

  ArrayAd<double, 1> od(nlayer, nparam);
  ArrayAd<double, 1> ssa(nlayer, nparam);
  
  for(int lay_idx = 0; lay_idx < nlayer; lay_idx++) {
    od(lay_idx) = taur(lay_idx) + taug(lay_idx);
    ssa(lay_idx) = taur(lay_idx) / od(lay_idx);
  }
  
  // No aerosols, and depolization factor = 0 
  // so simplified phase function moments:
  ArrayAd<double, 2> pf(3, nlayer, nparam);
  pf = 0.0;
  pf(0, all) = 1.0;
  pf(2, all) = 0.5;

  // Set up jacobians
  surface_params.jacobian() = 1.0;

  // Make copies of surface params as it will be modified
  // internally by the routines
  ArrayAd<double, 1> lid_surface_params;
  lid_surface_params = surface_params;
  ArrayAd<double, 1> ts_surface_params;
  ts_surface_params = surface_params;

  // Run lidort and 2stream to generate values for comparison
  lidort_driver.reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                surface_type, lid_surface_params,
                                                od, ssa, pf, refl_lid, jac_atm_lid, jac_surf_lid);
                                                
  twostream_driver.reflectance_and_jacobian_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, ts_surface_params,
                                                   od, ssa, pf, refl_ts, jac_atm_ts, jac_surf_ts);

  if(debug_output) {
    std::cerr << "refl_lid = " << refl_lid << std::endl
              << "refl_ts  = " << refl_ts << std::endl;
  }

  BOOST_CHECK_CLOSE(refl_lid, refl_ts, 7e-5);

  blitz::Array<double, 2> jac_atm_fd(nparam, nlayer);
  double refl_fd;

  for(int l_idx = 0; l_idx < nlayer; l_idx++) {
    for(int p_idx = 0; p_idx < nparam; p_idx++) {
      blitz::Array<double,1> taug_pert( taug.value().copy() );
      blitz::Array<double,1> taur_pert( taur.value().copy() );

      blitz::Array<double,1> od_pert( od.value().shape() );
      blitz::Array<double,1> ssa_pert( ssa.value().shape() );
      blitz::Array<double,2> pf_pert( pf.value().copy() );
   
      switch (p_idx) {
      case 0:
        taug_pert(l_idx) += pert_atm(p_idx);
        break;
      case 1:
        taur_pert(l_idx) += pert_atm(p_idx);
        break;
      }

      od_pert = taur_pert + taug_pert;
      ssa_pert = taur_pert / od_pert;

      refl_fd = twostream_driver.reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params.value().copy(),
                                                   od_pert, ssa_pert, pf_pert);

      jac_atm_fd(p_idx, l_idx) = (refl_fd - refl_ts) / pert_atm(p_idx);
    }
  }

  if(debug_output) {
    std::cerr << setprecision(8)
              << "jac_atm_lid = " << jac_atm_lid(rjac,rlay).transpose(1,0)
              << "jac_atm_fd = " << jac_atm_fd(rjac, rlay).transpose(1,0)
              << "jac_atm_ts = " << jac_atm_ts(rjac,rlay).transpose(1,0) << std::endl;
  }

  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_atm_lid(rjac,rlay), jac_atm_ts(rjac,rlay), 1e-6);
  BOOST_CHECK_MATRIX_CLOSE_TOL(jac_atm_ts(rjac,rlay), jac_atm_fd(rjac,rlay), 5e-4);

  if(blitz::any(surface_params.value() > 0.0)) {
    // Check surface jacobians against finite difference
    blitz::Array<double, 1> jac_surf_fd( jac_surf_ts.rows() );
    jac_surf_fd = 0.0;

    for(int p_idx = 0; p_idx < pert_surf.rows(); p_idx++) {
      blitz::Array<double,1> surface_params_pert( surface_params.rows() );
      surface_params_pert = surface_params.value();
      surface_params_pert(p_idx) += pert_surf(p_idx);

      refl_fd = twostream_driver.reflectance_calculate(heights, sza(0), zen(0), azm(0),
                                                   surface_type, surface_params_pert,
                                                   od.value(), ssa.value(), pf.value());

      jac_surf_fd(p_idx) = (refl_fd - refl_ts) / pert_surf(p_idx);

      // Adjust analytic jacobians to have same meaning as finite difference one
      jac_surf_lid(p_idx) *= lid_surface_params.jacobian()(p_idx, 0);
      jac_surf_ts(p_idx) *= ts_surface_params.jacobian()(p_idx, 0);
    }

    if(debug_output) {
      std::cerr << "jac_surf_lid = " << jac_surf_lid(rsurf) << std::endl
                << "jac_surf_fd = " << jac_surf_fd << std::endl
                << "jac_surf_ts = " << jac_surf_ts << std::endl;
    }

    BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_lid(rsurf), jac_surf_ts(rsurf), 1e-7);
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_ts, jac_surf_fd, 1e-7);
  }
}

BOOST_AUTO_TEST_CASE(lambertian)
{
  bool debug_output = false;

  int nlayer = 2;
  int nparam = 2;
  ArrayAd<double, 1> surface_params(1, 1); 
  ArrayAd<double, 1> taug(nlayer, nparam);
  ArrayAd<double, 1> taur(nlayer, nparam);

  // Use simple lambertian throughout
  int surface_type = LAMBERTIAN;
  
  // Check jacobian using finite derivatives
  // od, ssa
  blitz::Array<double, 1> pert_atm(2);

  // Perturbation value for surface fd checking
  Array<double, 1> pert_surf(1);

  ////////////////
  // Surface only
  surface_params.value() = 1.0;

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  if (debug_output) std::cerr << "Surface only" << std::endl << "----------------------------" << std::endl;  
  pert_atm = 1e-4, 1e-4;
  pert_surf = 1e-8;
  test_twostream(surface_type, surface_params, taug, taur, pert_atm, pert_surf, debug_output);

  ////////////////
  // Rayleigh only
  surface_params.value() = 1.0e-6;

  taur = 2.0e-2/nlayer;
  taug = 1.0e-6/nlayer;

  if (debug_output) std::cerr << "Rayleigh only" << std::endl << "----------------------------" << std::endl;  
  pert_atm = 1e-3, -1e-4;
  pert_surf = 1e-8;
  test_twostream(surface_type, surface_params, taug, taur, pert_atm, pert_surf, debug_output);

  ////////////////
  // Gas + Surface
  surface_params.value() = 1.0;

  taur = 1.0e-6/nlayer;
  taug = 1.0/nlayer;

  if (debug_output) std::cerr << "Gas + Surface" << std::endl << "----------------------------" << std::endl;  
  pert_atm = 1e-4, 1e-4;
  pert_surf = 1e-8;
  test_twostream(surface_type, surface_params, taug, taur, pert_atm, pert_surf, debug_output);
}

BOOST_AUTO_TEST_CASE(coxmunk)
{
  bool debug_output = false; 

  int nlayer = 2;
  int nparam = 2;
  ArrayAd<double, 1> surface_params(5, 4); 
  ArrayAd<double, 1> taug(nlayer, nparam);
  ArrayAd<double, 1> taur(nlayer, nparam);

  // Use simple lambertian throughout
  int surface_type = COXMUNK;
  
  // Check jacobian using finite derivatives
  // od, ssa
  blitz::Array<double, 1> pert_atm(2);
  pert_atm = 1e-4, 1e-4;
  
  // Surface only
  surface_params(0) = 1;
  surface_params(1) = 7;
  surface_params(2) = 1.334;
  surface_params(3) = 0.8;
  surface_params(4) = 0.0;

  // Surface perturbations
  blitz::Array<double, 1> pert_surf(surface_params.number_variable());
  pert_surf = 1e-4, 1e-4, 1e-4, 1e-4;
  //pert_surf(1) = sqrt(surface_params(1).value()*surface_params(1).value() + pert_surf(1)) - surface_params(1).value();

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  if (debug_output) std::cerr << "Coxmunk" << std::endl << "----------------------------" << std::endl;  
  test_twostream(surface_type, surface_params, taug, taur, pert_atm, pert_surf, debug_output);
}

BOOST_AUTO_TEST_CASE(brdf)
{
  bool debug_output = false;

  int nlayer = 2;
  int nparam = 2;
  ArrayAd<double, 1> surface_params(5, 1); 
  ArrayAd<double, 1> taug(nlayer, nparam);
  ArrayAd<double, 1> taur(nlayer, nparam);

  // Use simple lambertian throughout
  int surface_type = BREONVEG;
  
  // Check jacobian using finite derivatives
  // od, ssa
  blitz::Array<double, 1> pert_atm(2);

  // Perturbation value for surface fd checking
  Array<double, 1> pert_surf(5);
  pert_surf = 1e-8;

  ////////////////
  // Surface only
  surface_params.value()(0) = 1.0; // rahman kernel factor
  surface_params.value()(1) = 0.1; // hotspot parameter
  surface_params.value()(2) = 0.3; // asymmetry
  surface_params.value()(3) = 1.5; // anisotropy_parameter
  surface_params.value()(4) = 1.0; // breon kernel factor

  taur = 1.0e-6/nlayer;
  taug = 1.0e-6/nlayer;

  pert_atm = 1e-4, 1e-4;
  test_twostream(surface_type, surface_params, taug, taur, pert_atm, pert_surf, debug_output);
}

BOOST_AUTO_TEST_CASE(valgrind_problem)
{
  // This catchs a particular valgrind problem we encountered in the
  // full physics code. This originally came from wavenumber
  // 13050.47 in the first iteration of tccon_sounding_1. There are
  // other valgrind errors there, but the hope is fixing this one will
  // fix all of them.

  // Note this has been fixed, but we'll leave this test in place
  // for now.
  
  IfstreamCs captured_input(test_data_dir() + "expected/twostream/valgrind_problem");
  int nlayer, nparm, surface_type, do_fullquadrature;
  captured_input >> nlayer >> nparm >> surface_type 
                 >> do_fullquadrature;
  TwostreamRtDriver d(nlayer, surface_type, (do_fullquadrature == 1));
                      
  double refl;
  Array<double, 2> jac_atm;
  Array<double, 1> jac_surf;
  Array<double, 1> height;
  ArrayAd<double, 1> surface_param, od, ssa;
  ArrayAd<double, 2> pf;
  double sza, zen, azm;
  captured_input >> height >> sza >> azm >> zen >> surface_type
                 >> surface_param >> od >> ssa >> pf;
  d.reflectance_and_jacobian_calculate(height, sza, azm, zen, 
                                       surface_type, surface_param,
                                       od, ssa, pf, refl,
                                       jac_atm, jac_surf);
  // This will trigger an error when we run with valgrind. Note that
  // we *don't* actually see a NaN here, rather this conditional
  // triggers the unitialized value error
  for(int i = 0; i < od.jacobian().rows(); ++i)
    for(int j = 0; j < od.jacobian().cols(); ++j) {
      if(std::isnan(jac_atm(j,i)))
        std::cerr << "Nan at jac_atm(" << j << ", " << i << ")\n";
    }
  for(int i = 0; i < jac_surf.rows(); ++i)
    if(std::isnan(jac_surf(i)))
      std::cerr << "Nan at jac_surf(" << i << ")\n";
}
BOOST_AUTO_TEST_SUITE_END()
