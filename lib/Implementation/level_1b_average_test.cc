#include "unit_test_support.h"
#include "level_1b_acos.h"
#include "level_1b_average.h"
#include "gosat_noise_model.h"
#include "configuration_fixture.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(level_1b_average, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<HdfFile> hf(new HdfFile(test_data_dir() + "l1b.h5"));

  std::vector<boost::shared_ptr<Level1b> > l1b;
  boost::shared_ptr<AcosSoundingId> sid_p
    (new AcosSoundingId(*hf, "20090725020316", AcosSoundingId::P_SOUNDING));
  boost::shared_ptr<AcosSoundingId> sid_s
    (new AcosSoundingId(*hf, "20090725020316", AcosSoundingId::S_SOUNDING));
  boost::shared_ptr<NoiseModel> nm_p(new GosatNoiseModel(*hf, *sid_p, 
							 *config_instrument));
  boost::shared_ptr<Level1bAcos> l1b_p(new Level1bAcos(hf, sid_p));
  l1b_p->noise_model(nm_p);
  l1b.push_back(l1b_p);
  boost::shared_ptr<NoiseModel> nm_s(new GosatNoiseModel(*hf, *sid_s, 
							 *config_instrument));
  boost::shared_ptr<Level1bAcos> l1b_s(new Level1bAcos(hf, sid_s));
  l1b_s->noise_model(nm_s);
  l1b.push_back(l1b_s);
  Level1bAverage h(l1b);
  BOOST_CHECK_EQUAL(h.number_spectrometer(), 3);
  blitz::Array<double, 2> st_expect(3,4);
  st_expect = 
    1, -7.11679e-05, 0.00045423, -0.000161178,
    1, -0.000114352, 0.000830699, 0.000454567,
    1, -3.09348e-05, 0.000208143, 9.26852e-05;
  for(int i = 0; i < 3; ++i) {
    BOOST_CHECK_CLOSE(h.latitude(i).value, 67.30134, 1e-4);
    BOOST_CHECK_CLOSE(h.solar_zenith(i).value, 52.28602, 1e-4);
    BOOST_CHECK_CLOSE(h.solar_azimuth(i).value, 221.64485, 1e-4);
    BOOST_CHECK_CLOSE(h.altitude(i).value, 19.636714935302734, 1e-4);
    BOOST_CHECK_CLOSE(h.sounding_zenith(i).value, 15.066192626953125, 1e-4);
    BOOST_CHECK_CLOSE(h.sounding_azimuth(i).value, 283.14007568359375, 1e-4);
    BOOST_CHECK_MATRIX_CLOSE_TOL(h.stokes_coefficient(i), 
				 st_expect(i, Range::all()), 1e-5);
  }
  BOOST_CHECK_CLOSE(h.relative_velocity(0).value, -398.8081, 1e-4);
  BOOST_CHECK_CLOSE(h.time(0).pgs_time(), 5.226410045967183e8, 1e-4);
  BOOST_CHECK_EQUAL(h.sounding_id(), (int64_t) 20090725020316);
  BOOST_CHECK_EQUAL(h.exposure_index(), 152);

  BOOST_CHECK_EQUAL(h.radiance(0).data().extent(blitz::firstDim), 1805);
  BOOST_CHECK_CLOSE(h.radiance(0).data()(402 + 10), 5.2663784799733548e-07, 
		    1e-4);
  BOOST_CHECK_EQUAL(h.radiance(1).data().extent(blitz::firstDim), 3508);
  BOOST_CHECK_CLOSE(h.radiance(1).data()(2085 + 10), 4.7928787694218045e-07, 
		    1e-4);
  BOOST_CHECK_EQUAL(h.radiance(2).data().extent(blitz::firstDim), 2005);
  BOOST_CHECK_CLOSE(h.radiance(2).data()(300 + 10), 1.6280747416885788e-07, 
		    1e-4);
  IfstreamCs expected_data(test_data_dir() + "expected/level_1b_average/basic");
  Array<double, 1> rad_expected, uncert_expected;
  expected_data >> rad_expected;
  expected_data >> uncert_expected;
  BOOST_CHECK_MATRIX_CLOSE_TOL(h.radiance(0).data(), rad_expected, 1e-12);
  BOOST_CHECK_MATRIX_CLOSE_TOL(h.radiance(0).uncertainty(), uncert_expected, 
			       1e-14);
}

BOOST_AUTO_TEST_SUITE_END()
