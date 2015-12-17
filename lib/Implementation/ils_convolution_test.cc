#include "ils_convolution.h"
#include "hdf_file.h"
#include "dispersion_polynomial.h"
#include "ils_table.h"
#include "unit_test_support.h"
#include "ifstream_cs.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(ils_convolution, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Array<bool, 1> flag(2);
  flag = true, false;
  Array<double, 1> coeff(2);
  coeff = 1.28695614e+04, 1.99492886e-01;
  boost::shared_ptr<DispersionPolynomial>
    d(new DispersionPolynomial(coeff, flag, units::inv_cm, "Test band", 1805, 
			       true));
  HdfFile hf(test_data_dir() + "l2_fixed_level_static_input.h5");
  boost::shared_ptr<IlsTableLinear> ils_func(new IlsTableLinear(hf, 0, "A-Band", "o2"));
  IlsConvolution ils(d, ils_func);

  StateVector sv;
  sv.add_observer(ils);
  Array<double,1> x(2);
  x(0) = coeff(0);
  x(1) = 0;
  sv.update_state(x);

  std::vector<int> plist;
  // Arbitrary list of pixels.
  plist.push_back(403);
  plist.push_back(405);
  for(int i = 407; i <= 414; ++i)
    plist.push_back(i);

  IfstreamCs expected(test_data_dir() + "expected/ils_convolution/basic");
  Array<double, 1> wn_in, rad_hres_in, rad_out_expect;
  expected >> wn_in >> rad_hres_in >> rad_out_expect;

  BOOST_CHECK_MATRIX_CLOSE(ils.apply_ils(wn_in, rad_hres_in, plist),
			   rad_out_expect);

  Array<double, 2> jac_rad_fake(rad_hres_in.rows(), 2);
  jac_rad_fake = 0;
  jac_rad_fake(Range::all(), 1) = rad_hres_in;
  ArrayAd<double, 1> rad_hres_in2(rad_hres_in, jac_rad_fake);
  BOOST_CHECK_MATRIX_CLOSE(ils.apply_ils(wn_in, rad_hres_in2, plist).value(),
			   rad_out_expect);
  Array<double, 2> jac = ils.apply_ils(wn_in, rad_hres_in2, plist).jacobian();
  Array<double, 1> v0 = ils.apply_ils(wn_in, rad_hres_in2, plist).value();
  double epsilon = 1e-3;
  x(0) += epsilon;
  sv.update_state(x);
  Array<double, 1> v1 = ils.apply_ils(wn_in, rad_hres_in2, plist).value();
  Array<double, 2> jacd(jac.shape());
  jacd(Range::all(), 0) = (v1 - v0) / epsilon;
  x(0) -= epsilon;
  sv.update_state(x);
  rad_hres_in2.value() *= 1 + epsilon;
  Array<double, 1> v2 = ils.apply_ils(wn_in, rad_hres_in2, plist).value();
  jacd(Range::all(), 1) = (v2 - v0) / epsilon;
  BOOST_CHECK_MATRIX_CLOSE(jac, jacd);
}

BOOST_AUTO_TEST_SUITE_END()
