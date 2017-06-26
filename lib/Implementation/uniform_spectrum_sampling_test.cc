#include "uniform_spectrum_sampling.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(uniform_spectrum_sampling, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  UniformSpectrumSampling t(1.0,4.99,0.5);
  Array<double, 1> t_expect(9);
  t_expect = 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5;
  BOOST_CHECK_EQUAL(t.number_spectrometer(), 1);
  SpectralDomain dummy;		// Not used by UniformSpectrumSampling
  DoubleWithUnit dummy2;
  BOOST_CHECK_MATRIX_CLOSE(t.spectral_domain(0, dummy, dummy2).wavenumber(), 
			   t_expect);
}

BOOST_AUTO_TEST_CASE(threespect)
{
  UniformSpectrumSampling t(1.0,4.99,0.5,
			    1.0+1,4.99+1,0.5,
			    1.0+2,4.99+2,0.5
			    );
  Array<double, 1> t_expect(9);
  t_expect = 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5;
  BOOST_CHECK_EQUAL(t.number_spectrometer(), 3);
  SpectralDomain dummy;		// Not used by UniformSpectrumSampling
  DoubleWithUnit dummy2;
  BOOST_CHECK_MATRIX_CLOSE(t.spectral_domain(0, dummy, dummy2).wavenumber(), 
			   t_expect);
  BOOST_CHECK_MATRIX_CLOSE(t.spectral_domain(1, dummy, dummy2).wavenumber(), 
			   t_expect + 1);
  BOOST_CHECK_MATRIX_CLOSE(t.spectral_domain(2, dummy, dummy2).wavenumber(), 
			   t_expect + 2);
}

BOOST_AUTO_TEST_SUITE_END()
