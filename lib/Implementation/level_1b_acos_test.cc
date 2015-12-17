#include "unit_test_support.h"
#include "level_1b_acos.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(level_1b_acos, GlobalFixture)

BOOST_AUTO_TEST_CASE(p_type)
{
  boost::shared_ptr<HdfFile> hf(new HdfFile(test_data_dir() + "l1b.h5"));
  boost::shared_ptr<AcosSoundingId> sid
    (new AcosSoundingId(*hf, "20090725020316", AcosSoundingId::P_SOUNDING));
  Level1bAcos h(test_data_dir() + "l1b.h5", sid);
  Level1bAcos h2(hf, sid);
  BOOST_CHECK_EQUAL(h.number_spectrometer(), 3);
  BOOST_CHECK_EQUAL(h2.number_spectrometer(), 3);
  BOOST_CHECK(h.is_h_gain());
  BOOST_CHECK(h2.is_h_gain());
  BOOST_CHECK(!h.is_m_gain());
  BOOST_CHECK(!h2.is_m_gain());
  blitz::Array<double, 2> st_expect(3,4);
  st_expect = 
    1.00011, -0.987557, 0.0392983, -0.15301,
    1.00022, -0.975515, 0.0113057,  0.220646,
    1.00005, -0.982183, 0.0265311,  0.18633;
  for(int i = 0; i < 3; ++i) {
    BOOST_CHECK_CLOSE(h.latitude(i).value, 67.30134, 1e-4);
    BOOST_CHECK_CLOSE(h2.latitude(i).value, 67.30134, 1e-4);
    BOOST_CHECK_CLOSE(h.solar_zenith(i).value, 52.28602, 1e-4);
    BOOST_CHECK_CLOSE(h2.solar_zenith(i).value, 52.28602, 1e-4);
    BOOST_CHECK_CLOSE(h.solar_azimuth(i).value, 221.64485, 1e-4);
    BOOST_CHECK_CLOSE(h2.solar_azimuth(i).value, 221.64485, 1e-4);
    BOOST_CHECK_CLOSE(h.altitude(i).value, 19.636714935302734, 1e-4);
    BOOST_CHECK_CLOSE(h2.altitude(i).value, 19.636714935302734, 1e-4);
    BOOST_CHECK_CLOSE(h.sounding_zenith(i).value, 15.066192626953125, 1e-4);
    BOOST_CHECK_CLOSE(h2.sounding_zenith(i).value, 15.066192626953125, 1e-4);
    BOOST_CHECK_CLOSE(h.sounding_azimuth(i).value, 283.14007568359375, 1e-4);
    BOOST_CHECK_CLOSE(h2.sounding_azimuth(i).value, 283.14007568359375, 1e-4);
    BOOST_CHECK_CLOSE(h.land_fraction(i), 96.4072, 1e-4);
    BOOST_CHECK_CLOSE(h2.land_fraction(i), 96.4072, 1e-4);
    BOOST_CHECK_MATRIX_CLOSE_TOL(h.stokes_coefficient(i), 
				 st_expect(i, Range::all()), 1e-5);
    BOOST_CHECK_MATRIX_CLOSE_TOL(h2.stokes_coefficient(i), 
				 st_expect(i, Range::all()), 1e-5);
  }
  BOOST_CHECK_CLOSE(h.relative_velocity(0).value, -398.8081, 1e-4);
  BOOST_CHECK_CLOSE(h2.relative_velocity(0).value, -398.8081, 1e-4);
  BOOST_CHECK_CLOSE(h.time(0).pgs_time(), 5.226410045967183e8, 1e-4);
  BOOST_CHECK_CLOSE(h2.time(0).pgs_time(), 5.226410045967183e8, 1e-4);
  BOOST_CHECK_EQUAL(h.sounding_id(), (int64_t) 20090725020316);
  BOOST_CHECK_EQUAL(h2.sounding_id(), (int64_t) 20090725020316);
  BOOST_CHECK_EQUAL(h.exposure_index(), 152);
  BOOST_CHECK_EQUAL(h2.exposure_index(), 152);

  BOOST_CHECK_EQUAL(h.radiance(0).data().extent(blitz::firstDim), 1805);
  BOOST_CHECK_EQUAL(h2.radiance(0).data().extent(blitz::firstDim), 1805);
  BOOST_CHECK_CLOSE(h.radiance(0).data()(402 + 10), 5.0381521532472107e-07, 
		    1e-4);
  BOOST_CHECK_CLOSE(h2.radiance(0).data()(402 + 10), 5.0381521532472107e-07, 
		    1e-4);
  BOOST_CHECK_EQUAL(h.radiance(1).data().extent(blitz::firstDim), 3508);
  BOOST_CHECK_EQUAL(h2.radiance(1).data().extent(blitz::firstDim), 3508);
  BOOST_CHECK_CLOSE(h.radiance(1).data()(2085 + 10), 4.4335814664009376e-07, 
		    1e-4);
  BOOST_CHECK_CLOSE(h2.radiance(1).data()(2085 + 10), 4.4335814664009376e-07, 
		    1e-4);
  BOOST_CHECK_EQUAL(h.radiance(2).data().extent(blitz::firstDim), 2005);
  BOOST_CHECK_EQUAL(h2.radiance(2).data().extent(blitz::firstDim), 2005);
  BOOST_CHECK_CLOSE(h.radiance(2).data()(300 + 10), 1.4391470415375807e-07, 
		    1e-4);
  BOOST_CHECK_CLOSE(h2.radiance(2).data()(300 + 10), 1.4391470415375807e-07, 
		    1e-4);

  blitz::Array<double, 1> expt_coeff(2);
  expt_coeff = 12869.88457452, 0.19949289;
  BOOST_CHECK_MATRIX_CLOSE(h.spectral_coefficient(0).value, expt_coeff);
  expt_coeff =  5749.98346211, 0.19949289;
  BOOST_CHECK_MATRIX_CLOSE(h.spectral_coefficient(1).value, expt_coeff);
  expt_coeff =  4749.92562304, 0.19949289;
  BOOST_CHECK_MATRIX_CLOSE(h.spectral_coefficient(2).value, expt_coeff);
}

BOOST_AUTO_TEST_CASE(s_type)
{
  boost::shared_ptr<HdfFile> hf(new HdfFile(test_data_dir() + "l1b.h5"));
  boost::shared_ptr<AcosSoundingId> sid
    (new AcosSoundingId(*hf, "20090725020316", AcosSoundingId::S_SOUNDING));
  Level1bAcos h(test_data_dir() + "l1b.h5", sid);
  Level1bAcos h2(hf, sid);
  BOOST_CHECK_EQUAL(h.number_spectrometer(), 3);
  blitz::Array<double, 2> st_expect(3,4);
  st_expect = 
    0.999888, 0.987415, -0.0383899,   0.152688,
    0.99978,  0.975286, -0.00964431, -0.219737,
    0.999947, 0.982121, -0.0261149,  -0.186145;
  for(int i = 0; i < 3; ++i) {
    BOOST_CHECK_CLOSE(h.latitude(i).value, 67.30134, 1e-4);
    BOOST_CHECK_CLOSE(h2.latitude(i).value, 67.30134, 1e-4);
    BOOST_CHECK_CLOSE(h.solar_zenith(i).value, 52.28602, 1e-4);
    BOOST_CHECK_CLOSE(h2.solar_zenith(i).value, 52.28602, 1e-4);
    BOOST_CHECK_CLOSE(h.solar_azimuth(i).value, 221.64485, 1e-4);
    BOOST_CHECK_CLOSE(h2.solar_azimuth(i).value, 221.64485, 1e-4);
    BOOST_CHECK_CLOSE(h.altitude(i).value, 19.636714935302734, 1e-4);
    BOOST_CHECK_CLOSE(h2.altitude(i).value, 19.636714935302734, 1e-4);
    BOOST_CHECK_CLOSE(h.sounding_zenith(i).value, 15.066192626953125, 1e-4);
    BOOST_CHECK_CLOSE(h2.sounding_zenith(i).value, 15.066192626953125, 1e-4);
    BOOST_CHECK_CLOSE(h.sounding_azimuth(i).value, 283.14007568359375, 1e-4);
    BOOST_CHECK_CLOSE(h2.sounding_azimuth(i).value, 283.14007568359375, 1e-4);
    BOOST_CHECK_CLOSE(h.land_fraction(i), 96.4072, 1e-4);
    BOOST_CHECK_CLOSE(h2.land_fraction(i), 96.4072, 1e-4);
    BOOST_CHECK_MATRIX_CLOSE_TOL(h.stokes_coefficient(i), 
				 st_expect(i, Range::all()), 1e-5);
    BOOST_CHECK_MATRIX_CLOSE_TOL(h2.stokes_coefficient(i), 
				 st_expect(i, Range::all()), 1e-5);
  }
  BOOST_CHECK_CLOSE(h.relative_velocity(0).value, -398.8081, 1e-4);
  BOOST_CHECK_CLOSE(h2.relative_velocity(0).value, -398.8081, 1e-4);
  BOOST_CHECK_CLOSE(h.time(0).pgs_time(), 5.226410045967183e8, 1e-4);
  BOOST_CHECK_CLOSE(h2.time(0).pgs_time(), 5.226410045967183e8, 1e-4);
  BOOST_CHECK_EQUAL(h.sounding_id(), (int64_t) 20090725020316);
  BOOST_CHECK_EQUAL(h2.sounding_id(), (int64_t) 20090725020316);
  BOOST_CHECK_EQUAL(h.exposure_index(), 152);
  BOOST_CHECK_EQUAL(h2.exposure_index(), 152);

  BOOST_CHECK_EQUAL(h.radiance(0).data().extent(blitz::firstDim), 1805);
  BOOST_CHECK_EQUAL(h2.radiance(0).data().extent(blitz::firstDim), 1805);
  BOOST_CHECK_CLOSE(h.radiance(0).data()(402 + 10), 5.494604806699499e-07, 
		    1e-4);
  BOOST_CHECK_CLOSE(h2.radiance(0).data()(402 + 10), 5.494604806699499e-07, 
		    1e-4);
  BOOST_CHECK_EQUAL(h.radiance(1).data().extent(blitz::firstDim), 3508);
  BOOST_CHECK_EQUAL(h2.radiance(1).data().extent(blitz::firstDim), 3508);
  BOOST_CHECK_CLOSE(h.radiance(1).data()(2085 + 10), 5.1521760724426713e-07, 
		    1e-4);
  BOOST_CHECK_CLOSE(h2.radiance(1).data()(2085 + 10), 5.1521760724426713e-07, 
		    1e-4);
  BOOST_CHECK_EQUAL(h.radiance(2).data().extent(blitz::firstDim), 2005);
  BOOST_CHECK_EQUAL(h2.radiance(2).data().extent(blitz::firstDim), 2005);
  BOOST_CHECK_CLOSE(h.radiance(2).data()(300 + 10), 1.8170024418395769e-07, 
		    1e-4);
  BOOST_CHECK_CLOSE(h2.radiance(2).data()(300 + 10), 1.8170024418395769e-07, 
		    1e-4);

  blitz::Array<double, 1> expt_coeff(2);
  expt_coeff = 12869.88457452, 0.19949289;
  BOOST_CHECK_MATRIX_CLOSE(h.spectral_coefficient(0).value, expt_coeff);
  expt_coeff =  5749.98346211, 0.19949289;
  BOOST_CHECK_MATRIX_CLOSE(h.spectral_coefficient(1).value, expt_coeff);
  expt_coeff =  4749.92562304, 0.19949289;
  BOOST_CHECK_MATRIX_CLOSE(h.spectral_coefficient(2).value, expt_coeff);
}


BOOST_AUTO_TEST_CASE(bad_data)
{
  boost::shared_ptr<HdfFile> hf(new HdfFile(test_data_dir() + "l1b.h5"));
  boost::shared_ptr<AcosSoundingId> sid
    (new AcosSoundingId(*hf, "20090725020316", AcosSoundingId::P_SOUNDING));
  Level1bAcos h(hf, sid);
  BOOST_CHECK_THROW(h.latitude(-1), Exception);
  BOOST_CHECK_THROW(h.latitude(3), Exception);
  BOOST_CHECK_THROW(h.solar_zenith(-1), Exception);
  BOOST_CHECK_THROW(h.solar_zenith(3), Exception);
  BOOST_CHECK_THROW(h.solar_azimuth(-1), Exception);
  BOOST_CHECK_THROW(h.solar_azimuth(3), Exception);
  BOOST_CHECK_THROW(h.altitude(-1), Exception);
  BOOST_CHECK_THROW(h.altitude(3), Exception);
  BOOST_CHECK_THROW(h.land_fraction(-1), Exception);
  BOOST_CHECK_THROW(h.land_fraction(3), Exception);

}

BOOST_AUTO_TEST_SUITE_END()
