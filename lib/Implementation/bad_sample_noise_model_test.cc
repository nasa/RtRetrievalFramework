#include "unit_test_support.h"
#include "bad_sample_noise_model.h"
#include "oco_noise_model.h"
#include "oco_sounding_id.h"
#include "level_1b_oco.h"
#include "fp_exception.h"
#include <iostream>

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(bad_sample_noise_model, GlobalFixture)
BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<HdfFile> h(new HdfFile(test_data_dir() + "/oco2_L1bScND_80008a_111017225030d_spliced.h5"));
  boost::shared_ptr<OcoSoundingId> sid(new OcoSoundingId(*h, "2010090900133834"));
  Level1bOco l1b = Level1bOco(h, sid);
  boost::shared_ptr<NoiseModel> noise_model(new OcoNoiseModel(*h, *sid));

  Array<double,1> uncertainty_calc = 
    noise_model->uncertainty(0, l1b.radiance(0).data());
  Array<bool, 2> bad_sample(3, uncertainty_calc.rows());
  bad_sample = false;
  bad_sample(0, 10) = true;
  BadSampleNoiseModel bnm(noise_model, bad_sample, 1e20);
  uncertainty_calc(10) = 1e20;
  BOOST_CHECK_MATRIX_CLOSE_TOL(bnm.uncertainty(0, l1b.radiance(0).data()),
			       uncertainty_calc, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

