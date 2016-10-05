#include "unit_test_support.h"
#include "gosat_noise_model.h"
#include "configuration_fixture.h"
#include "acos_sounding_id.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(gosat_noise_model, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(p_type_hdf)
{
  HdfFile hdf_file(test_data_dir() + "l1b.h5");
  AcosSoundingId snd_id(hdf_file, "20090725020316", AcosSoundingId::P_SOUNDING);
  GosatNoiseModel noise_model(hdf_file, snd_id, *config_instrument);

  IfstreamCs expected_data(test_data_dir() + "expected/gosat_noise_model/p_radiance_uncertainty");
  Array<double, 1> exp_rad, exp_uncert;
  expected_data >> exp_rad;
  expected_data >> exp_uncert;

  Array<double, 1> calc_uncert = noise_model.uncertainty(0, exp_rad);
  
  BOOST_CHECK_MATRIX_CLOSE_TOL(exp_uncert, calc_uncert, 1e-14);
}

BOOST_AUTO_TEST_CASE(s_type_hdf)
{
  HdfFile hdf_file(test_data_dir() + "l1b.h5");
  AcosSoundingId snd_id(hdf_file, "20090725020316", AcosSoundingId::S_SOUNDING);
  GosatNoiseModel noise_model(hdf_file, snd_id, *config_instrument);

  IfstreamCs expected_data(test_data_dir() + "expected/gosat_noise_model/s_radiance_uncertainty");
  Array<double, 1> exp_rad, exp_uncert;
  expected_data >> exp_rad;
  expected_data >> exp_uncert;

  Array<double, 1> calc_uncert = noise_model.uncertainty(0, exp_rad);
  
  BOOST_CHECK_MATRIX_CLOSE_TOL(exp_uncert, calc_uncert, 1e-14);
}

BOOST_AUTO_TEST_CASE(ascii_read)
{
  HeritageFile noise_ascii(test_data_dir() + "in/noise_cnv_avg_H.dat");
  GosatNoiseModel noise_model(noise_ascii);

  IfstreamCs expected_data(test_data_dir() + "expected/gosat_noise_model/avg_radiance_uncertainty");
  Array<double, 1> exp_rad, exp_uncert;
  expected_data >> exp_rad;
  expected_data >> exp_uncert;

  Array<double, 1> calc_uncert = noise_model.uncertainty(0, exp_rad);

  BOOST_CHECK_MATRIX_CLOSE_TOL(exp_uncert, calc_uncert, 1e-14);
}

BOOST_AUTO_TEST_CASE(emp_noise_hdf)
{
  HdfFile hdf_file(test_data_dir() + "l1b.h5");
  HeritageFile emp_noise_file(test_data_dir() + "in/emp_noise_coeff.dat");
  AcosSoundingId snd_id(hdf_file, "20090725020316", AcosSoundingId::P_SOUNDING);
  GosatNoiseModel noise_model(hdf_file, snd_id, *config_instrument, 
			      emp_noise_file);

  IfstreamCs expected_data(test_data_dir() + "expected/gosat_noise_model/emp_noise_radiance_uncertainty");
  Array<double, 1> exp_rad, exp_uncert;
  expected_data >> exp_rad;
  expected_data >> exp_uncert;

  Array<double, 1> calc_uncert = noise_model.uncertainty(0, exp_rad);

  BOOST_CHECK_MATRIX_CLOSE_TOL(exp_uncert, calc_uncert, 1e-14);
}

BOOST_AUTO_TEST_CASE(emp_noise_hdf2)
{
  HdfFile hdf_file(test_data_dir() + "l1b.h5");
  HdfFile hdf_static_input(test_data_dir() + "l2_fixed_level_static_input.h5");
  AcosSoundingId snd_id(hdf_file, "20090725020316", AcosSoundingId::P_SOUNDING);
  GosatNoiseModel noise_model(hdf_file, snd_id, *config_instrument, 
			      hdf_static_input, "Instrument");

  IfstreamCs expected_data(test_data_dir() + "expected/gosat_noise_model/emp_noise_radiance_uncertainty");
  Array<double, 1> exp_rad, exp_uncert;
  expected_data >> exp_rad;
  expected_data >> exp_uncert;

  Array<double, 1> calc_uncert = noise_model.uncertainty(0, exp_rad);
  
  BOOST_CHECK_MATRIX_CLOSE_TOL(exp_uncert, calc_uncert, 1e-14);
}


BOOST_AUTO_TEST_SUITE_END()
