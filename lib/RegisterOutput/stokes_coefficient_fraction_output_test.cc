#include "unit_test_support.h"
#include "stokes_coefficient_fraction_output.h"
#include "output_hdf.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(stokes_coefficient_fraction_output, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Array<double, 2> coeff(3, 4);
  coeff = 
    1, 2, 3, 4,
    5, 6, 7, 8,
    9, 10, 11, 12;
  Array<double, 1> f(3);
  Array<bool, 1> used(3);
  f = 0.1, 0.2, 0.3;
  used = true, true, true;
  boost::shared_ptr<StokesCoefficientFraction> s
    (new StokesCoefficientFraction(coeff, f, used));
  
  StokesCoefficientFractionOutput po(s, 0, "o2");
  boost::shared_ptr<OutputHdf> out(new OutputHdf("stokes_coefficient_fraction_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("stokes_coefficient_fraction_output.h5");
  po.register_output_apriori(out);
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in gosat instrument unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


