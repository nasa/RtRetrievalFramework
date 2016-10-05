#include "spectral_domain.h"
#include "fp_exception.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(spectral_domain, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  Array<double, 1> wavenumber(5);
  wavenumber = 1000, 1000.5, 1001, 1001.5, 1002;
  Array<double, 1> wavelength(5);
  wavelength = 1e4 / wavenumber;
  SpectralDomain wn(wavenumber);
  SpectralDomain wl(wavelength, units::micron);
  BOOST_CHECK_MATRIX_CLOSE(wn.data(), wavenumber);
  BOOST_CHECK_EQUAL(wn.units().name(), "cm^-1");
  BOOST_CHECK_EQUAL(wn.type_preference(), SpectralDomain::PREFER_WAVENUMBER);
  BOOST_CHECK_MATRIX_CLOSE(wn.wavenumber(), wavenumber);
  BOOST_CHECK_MATRIX_CLOSE(wn.wavelength(), wavelength);
  BOOST_CHECK_MATRIX_CLOSE(wn.wavenumber(1 / units::m), wavenumber * 100.0);
  BOOST_CHECK_MATRIX_CLOSE(wn.wavelength(units::nm), wavelength * 1e3);
  BOOST_CHECK_MATRIX_CLOSE(wl.data(), wavelength);
  BOOST_CHECK_EQUAL(wl.units().name(), "micron");
  BOOST_CHECK_EQUAL(wl.type_preference(), SpectralDomain::PREFER_WAVELENGTH);
  BOOST_CHECK_MATRIX_CLOSE(wl.wavenumber(), wavenumber);
  BOOST_CHECK_MATRIX_CLOSE(wl.wavelength(), wavelength);
  BOOST_CHECK_MATRIX_CLOSE(wl.wavenumber(1 / units::m), wavenumber * 100.0);
  BOOST_CHECK_MATRIX_CLOSE(wl.wavelength(units::nm), wavelength * 1e3);
}

BOOST_AUTO_TEST_SUITE_END()

