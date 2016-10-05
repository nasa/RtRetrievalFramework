#include "unit_test_support.h"
#include "dispersion_polynomial_output.h"
#include "output_hdf.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(dispersion_polynomial_output, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Array<bool, 1> flag(2);
  flag = true, true;
  Array<double, 1> coeff(2);
  coeff = 1.28695614e+04, 1.99492886e-01;
  boost::shared_ptr<DispersionPolynomial> 
    d(new DispersionPolynomial(coeff, flag, units::inv_cm, "Test band", 1805, 
			       true));
  
  DispersionPolynomialOutput po(d, "o2");
  boost::shared_ptr<OutputHdf> out(new OutputHdf("dispersion_polynomial_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("dispersion_polynomial_output.h5");
  po.register_output_apriori(out);
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in gosat instrument unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


