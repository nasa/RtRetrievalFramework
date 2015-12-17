#include "unit_test_support.h"
#include "empirical_orthogonal_function_output.h"
#include "ils_instrument.h"
#include "ils_convolution.h"
#include "hdf_file.h"
#include "dispersion_polynomial.h"
#include "ils_table.h"
#include "output_hdf.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(empirical_orthogonal_function_output, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > > corr(3);

  HdfFile hf(test_data_dir() + "l2_fixed_level_static_input.h5");
  Array<bool, 1> flag(2);
  flag = true, false;
  Array<double, 1> coeff(2);
  coeff = 1.28695614e+04, 1.99492886e-01;
  boost::shared_ptr<IlsTableLinear> 
    ils_tab(new IlsTableLinear(hf, 0, "A-Band", "o2"));
  boost::shared_ptr<DispersionPolynomial>
    d(new DispersionPolynomial(coeff, flag, units::inv_cm, 
			       ils_tab->band_name(), 1805, true));
  std::vector<boost::shared_ptr<Ils> > ils;
  ils.push_back(boost::shared_ptr<Ils>(new IlsConvolution(d, ils_tab)));
  boost::shared_ptr<EmpiricalOrthogonalFunction> 
    zoff(new EmpiricalOrthogonalFunction(1.0, true, *d, hf, 0, 0, 1, ils_tab->band_name()));

  EmpiricalOrthogonalFunctionOutput po(zoff, "o2");
  boost::shared_ptr<OutputHdf> out(new OutputHdf("empirical_orthogonal_function_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("empirical_orthogonal_function_output.h5");
  po.register_output_apriori(out);
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in gosat instrument unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


