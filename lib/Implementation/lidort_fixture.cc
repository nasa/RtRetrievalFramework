#include "lidort_fixture.h"
#include "stokes_coefficient_constant.h"
#include "spurr_brdf_types.h"

using namespace FullPhysics;
using namespace blitz;

LidortDriverCommonFixture::LidortDriverCommonFixture() : sza(3), zen(3), azm(3)
{
  // Set viewing geometry
  sza = 0.001;
  zen = 0.0;
  azm = 0.0;
  pure_nadir = false;

  wn_arr.resize(2);
  wn_arr = 12929.94, 12930.30; 
}

LidortDriverLambertianFixture::LidortDriverLambertianFixture() : LidortDriverCommonFixture()
{
  int nstreams = 4;
  int nmoms = 3;
  bool do_multiple_scattering_only = false;

  lidort_driver.reset(new LidortRtDriver(nstreams, nmoms, do_multiple_scattering_only, LAMBERTIAN, zen, pure_nadir));  
}

LidortDriverCoxmunkFixture::LidortDriverCoxmunkFixture() : LidortDriverCommonFixture()
{
  int nstreams = 4;
  int nmoms = 2*nstreams;
  bool do_multiple_scattering_only = false;

  lidort_driver.reset(new LidortRtDriver(nstreams, nmoms, do_multiple_scattering_only, COXMUNK, zen, pure_nadir));
}

LidortRtCommonFixture::LidortRtCommonFixture() : sza(3), zen(3), azm(3)
{
  // Set viewing geometry
  sza = 74.128288268999995;
  zen = 1.0e-6;
  azm = 12.504928589999992;
  pure_nadir = false;

  wn_arr.resize(2);
  wn_arr = 12929.94, 12930.30; 

  blitz::Array<double, 2> stokes_coef_v(3, 3);
  stokes_coef_v = 
    1,0,0,
    1,0,0,
    1,0,0;
  stokes_coefs.reset(new StokesCoefficientConstant(stokes_coef_v));
}

LidortLambertianFixture::LidortLambertianFixture(const std::string& Config_file) : 
  LidortRtCommonFixture(), ConfigurationFixture(Config_file)
{
  int nstreams = 4;
  int nmoms = 2*nstreams;
  bool do_multiple_scattering_only = false;
  
  lidort_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			       nstreams, nmoms, do_multiple_scattering_only));
  
  pp_expected_filename = test_data_dir() + "expected/lidort_driver/lambertian_plane_parallel";
  ps_expected_filename = test_data_dir() + "expected/lidort_driver/lambertian_pseudo_spherical";
  pp_and_ss_expected_filename = test_data_dir() + "expected/lidort_driver/lambertian_pp_plus_sscorrection";
}

LidortCoxmunkFixture::LidortCoxmunkFixture(const std::string& Config_file) : 
  LidortRtCommonFixture(), ConfigurationCoxmunkFixture(Config_file)
{
  int nstreams = 4;
  int nmoms = 2*nstreams;
  bool do_multiple_scattering_only = false;

  lidort_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			       nstreams, nmoms, do_multiple_scattering_only));

  pp_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk_plane_parallel";
  ps_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk_pseudo_spherical";
  pp_and_ss_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk_pp_plus_sscorrection";

}

LidortCoxmunkPlusLambertianFixture::LidortCoxmunkPlusLambertianFixture() : LidortRtCommonFixture()
{
  int nstreams = 4;
  int nmoms = 2*nstreams;
  bool do_multiple_scattering_only = false;

  lidort_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			       nstreams, nmoms, do_multiple_scattering_only));

  pp_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk+lamb_plane_parallel";
  ps_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk+lamb_pseudo_spherical";
  pp_and_ss_expected_filename = test_data_dir() + "expected/lidort_driver/coxmunk+lamb_pp_plus_sscorrection";

}

/*******************************************************************/

LidortLowHighCommon::LidortLowHighCommon() 
: sza(3), zen(3), azm(3)
{
  // Set viewing geometry
  sza = 74.128288268999995;
  zen = 0;
  azm = 12.504928589999992;
  pure_nadir = true;

  blitz::Array<double, 2> stokes_coef_v(3, 3);
  stokes_coef_v = 
    1,0,0,
    1,0,0,
    1,0,0;
  stokes_coefs.reset(new StokesCoefficientConstant(stokes_coef_v));

  nstreams_low = 2;
  nmoms_low = 2*nstreams_low;

  nstreams_high = 8;
  nmoms_high = Lidort_Pars::instance().maxmoments_input;

}

LidortLowHighLambertianFixture::LidortLowHighLambertianFixture() : LidortLowHighCommon()
{
  bool do_multiple_scattering_only = true;

  low_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			    nstreams_low, nmoms_low, do_multiple_scattering_only));

  high_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			     nstreams_high, nmoms_high, do_multiple_scattering_only));

}

LidortLowHighFixtureNoPolarization::LidortLowHighFixtureNoPolarization() : LidortLowHighCommon()
{
  bool do_multiple_scattering_only = false;

  low_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			    nstreams_low, nmoms_low, do_multiple_scattering_only));

  high_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			     nstreams_high, nmoms_high, do_multiple_scattering_only));
}

LidortLowHighCoxmunkFixture::LidortLowHighCoxmunkFixture() : LidortLowHighCommon()
{
  bool do_multiple_scattering_only = true;

  low_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			    nstreams_low, nmoms_low, do_multiple_scattering_only));

  high_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			     nstreams_high, nmoms_high, do_multiple_scattering_only));
}

LidortLowHighCoxmunkPlusLambertianFixture::LidortLowHighCoxmunkPlusLambertianFixture() : LidortLowHighCommon()
{
  bool do_multiple_scattering_only = true;

  low_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			    nstreams_low, nmoms_low, do_multiple_scattering_only));

  high_rt.reset(new LidortRt(config_atmosphere, stokes_coefs, sza, zen, azm, pure_nadir,
			     nstreams_high, nmoms_high, do_multiple_scattering_only));
}
