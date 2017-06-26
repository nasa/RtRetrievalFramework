#include "unit_test_support.h"
#include "oco_noise_model.h"
#include "oco_sounding_id.h"
#include "level_1b_oco.h"
#include "fp_exception.h"
#include <iostream>

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(oco_noise_model, GlobalFixture)

BOOST_AUTO_TEST_CASE(old_max_ms)
{
  boost::shared_ptr<HdfFile> h(new HdfFile(test_data_dir() + "/oco2_L1bScND_80008a_111017225030d_spliced.h5"));
  boost::shared_ptr<OcoSoundingId> sid(new OcoSoundingId(*h, "2010090900133834"));
  Level1bOco l1b = Level1bOco(h, sid);

  Array<double, 1> old_max_ms(3);
  old_max_ms = 1.4e21, 4.9e20, 1.7e20;
  OcoNoiseModel noise_model = OcoNoiseModel(*h, *sid, old_max_ms);

  Array<double, 2> expt_p_beg_end(3, 2);
  expt_p_beg_end = 
    0.010034451228113186,  0.0085836252818802059,
    0.0076201586316476388, 0.0083756448403055593,
    0.0086108141077417422, 0.010608880196296272;

  Array<double, 2> expt_b_beg_end(3, 2);
  expt_b_beg_end = 
    0.017036250081369633,  0.020973428336803795,
    0.0093927604913324053, 0.010068064453519163,
    0.029241598560507613,  0.03859064426935771;

 for(int s_idx = 0; s_idx < 3; s_idx++) {
    BOOST_CHECK_CLOSE(expt_p_beg_end(s_idx, 0), noise_model.coef_photon(s_idx)(0), 1e-8);
    BOOST_CHECK_CLOSE(expt_p_beg_end(s_idx, 1), noise_model.coef_photon(s_idx)(1015), 1e-8);
    BOOST_CHECK_CLOSE(expt_b_beg_end(s_idx, 0), noise_model.coef_background(s_idx)(0), 1e-8);
    BOOST_CHECK_CLOSE(expt_b_beg_end(s_idx, 1), noise_model.coef_background(s_idx)(1015), 1e-8);
  }

  Array<double,1> uncertainty_calc = noise_model.uncertainty(0, l1b.radiance(0).data()); 

  IfstreamCs expt_file(test_data_dir() + "expected/oco_noise_model/old_maxms");
  Array<double,1> uncertainty_expt;
  expt_file >> uncertainty_expt;
  // These are big numbers (1e20) so need a higher tolerance value
  BOOST_CHECK_MATRIX_CLOSE_TOL(uncertainty_expt, uncertainty_calc, 1e3);
}

BOOST_AUTO_TEST_CASE(hdf_read_compute)
{
  boost::shared_ptr<HdfFile> h(new HdfFile(test_data_dir() + "/oco2_L1bScND_80008a_111017225030d_spliced.h5"));
  boost::shared_ptr<OcoSoundingId> sid(new OcoSoundingId(*h, "2010090900133834"));
  Level1bOco l1b = Level1bOco(h, sid);
  OcoNoiseModel noise_model = OcoNoiseModel(*h, *sid);

  Array<double, 2> expt_p_beg_end(3, 2);
  expt_p_beg_end = 
    0.010034451228113186,  0.0085836252818802059,
    0.0076201586316476388, 0.0083756448403055593,
    0.0086108141077417422, 0.010608880196296272;

  Array<double, 2> expt_b_beg_end(3, 2);
  expt_b_beg_end = 
    0.017036250081369633,  0.020973428336803795,
    0.0093927604913324053, 0.010068064453519163,
    0.029241598560507613,  0.03859064426935771;

 for(int s_idx = 0; s_idx < 3; s_idx++) {
    BOOST_CHECK_CLOSE(expt_p_beg_end(s_idx, 0), noise_model.coef_photon(s_idx)(0), 1e-8);
    BOOST_CHECK_CLOSE(expt_p_beg_end(s_idx, 1), noise_model.coef_photon(s_idx)(1015), 1e-8);
    BOOST_CHECK_CLOSE(expt_b_beg_end(s_idx, 0), noise_model.coef_background(s_idx)(0), 1e-8);
    BOOST_CHECK_CLOSE(expt_b_beg_end(s_idx, 1), noise_model.coef_background(s_idx)(1015), 1e-8);
  }

  Array<double,1> uncertainty_calc = noise_model.uncertainty(0, l1b.radiance(0).data()); 

  IfstreamCs expt_file(test_data_dir() + "expected/oco_noise_model/hdf_file");
  Array<double,1> uncertainty_expt;
  expt_file >> uncertainty_expt;
  // These are big numbers (1e20) so need a higher tolerance value
  BOOST_CHECK_MATRIX_CLOSE_TOL(uncertainty_expt, uncertainty_calc, 1e3);
}

BOOST_AUTO_TEST_SUITE_END()

