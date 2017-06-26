#include "precomputed_noise_model.h"
#include "unit_test_support.h"
#include <iostream>

#include <boost/assign/std/vector.hpp> // for vector 'operator+=()'

using namespace FullPhysics;
using namespace blitz;
using namespace boost::assign; // bring 'operator+=()' into scope

BOOST_FIXTURE_TEST_SUITE(precomputed_noise_model, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  PrecomputedNoiseModel noise(test_data_dir() + "in/l1b/spec/spectra.dat");

  std::vector< Array<double,1> > empty_rad;
  empty_rad += Array<double,1>(1805), Array<double,1>(3508), Array<double,1>(2005);

  std::vector< Array<double,1> > exp_uncerts(3);
  IfstreamCs exp_file(test_data_dir() + "expected/precomputed_noise_model/uncertainty");
  for(unsigned int spec_idx = 0; spec_idx < exp_uncerts.size(); spec_idx++) {
    exp_file >> exp_uncerts[spec_idx]; 
  }

  for (unsigned int spec_idx = 0; spec_idx < empty_rad.size(); spec_idx++) {
    Array<double, 1> read_uncert(noise.uncertainty(spec_idx, empty_rad[spec_idx]));
    BOOST_CHECK_MATRIX_CLOSE_TOL(exp_uncerts[spec_idx], read_uncert, 1e-14);
  }

}

BOOST_AUTO_TEST_SUITE_END()
