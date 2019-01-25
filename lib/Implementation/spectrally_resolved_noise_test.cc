#include "unit_test_support.h"
#include "spectrally_resolved_noise.h"

#include "gosat_noise_model.h"
#include <boost/assign/std/vector.hpp> // for vector 'operator+=()'

using namespace FullPhysics;
using namespace blitz;
using namespace boost::assign; // bring 'operator+=()' into scope

BOOST_FIXTURE_TEST_SUITE(spectrally_resolved_noise, GlobalFixture)

BOOST_AUTO_TEST_CASE(base)
{
  HeritageFile noise_ascii(test_data_dir() + "in/noise_cnv_avg_H.dat");
  boost::shared_ptr<GosatNoiseModel> base_noise(new GosatNoiseModel(noise_ascii));
  SpectrallyResolvedNoise spec_noise(base_noise);

  std::vector< Array<double,1> > empty_rad;
  empty_rad += Array<double,1>(1805), Array<double,1>(3508), Array<double,1>(2005);

  // Check that uncertainty is passed through if no coefficents are set
  for (unsigned int spec_idx = 0; spec_idx < empty_rad.size(); spec_idx++) {
    Array<double, 1> base_uncert(base_noise->uncertainty(spec_idx, empty_rad[spec_idx]));
    Array<double, 1> spec_uncert(spec_noise.uncertainty(spec_idx, empty_rad[spec_idx]));
    BOOST_CHECK_MATRIX_CLOSE_TOL(base_uncert, spec_uncert, 1e-14);
  }

  for (unsigned int spec_idx = 0; spec_idx < empty_rad.size(); spec_idx++) {
    double scaling = (spec_idx + 1) / 100.0;
    Array<double, 1> base_uncert(base_noise->uncertainty(spec_idx, empty_rad[spec_idx]));
    base_uncert = base_uncert * scaling;

    // Check that scaling works when applied across the whole band
    Array<double, 1> coeffs(empty_rad[spec_idx].rows());
    coeffs = scaling;
    spec_noise.set_full_noise_scaling(spec_idx, coeffs);

    Array<double, 1> spec_uncert(spec_noise.uncertainty(spec_idx, empty_rad[spec_idx]));
    BOOST_CHECK_MATRIX_CLOSE_TOL(base_uncert, spec_uncert, 1e-14);

    // Check that a single scaling value works
    spec_noise.set_single_noise_scaling(spec_idx, scaling);

    spec_uncert.reference(spec_noise.uncertainty(spec_idx, empty_rad[spec_idx]));
    BOOST_CHECK_MATRIX_CLOSE_TOL(base_uncert, spec_uncert, 1e-14);

  }

}

BOOST_AUTO_TEST_SUITE_END()
