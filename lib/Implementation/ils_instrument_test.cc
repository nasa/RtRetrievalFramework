#include "ils_instrument.h"
#include "ils_convolution.h"
#include "hdf_file.h"
#include "dispersion_polynomial.h"
#include "ils_table.h"
#include "unit_test_support.h"
#include "configuration_fixture.h"
#include "accumulated_timer.h"
#include <iostream>
using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(ils_instrument, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  HdfFile hf(test_data_dir() + "l2_fixed_level_static_input.h5");
  Array<bool, 1> flag(2);
  flag = true, false;
  Array<double, 1> coeff(2);
  coeff = 1.28695614e+04, 1.99492886e-01;
  boost::shared_ptr<IlsTableLinear> ils_tab(new IlsTableLinear(hf, 0, "A-Band", "o2"));
  boost::shared_ptr<DispersionPolynomial>
    d(new DispersionPolynomial(coeff, flag, units::inv_cm, ils_tab->band_name(),
			       1805, true));
  std::vector<boost::shared_ptr<Ils> > ils;
  ils.push_back(boost::shared_ptr<Ils>(new IlsConvolution(d, ils_tab)));
  
  coeff = 5.74982835e+03, 1.99492886e-01;
  ils_tab.reset(new IlsTableLinear(hf, 1, "WC-Band", "weak_co2"));
  d.reset(new DispersionPolynomial(coeff, flag, units::inv_cm, 
				   ils_tab->band_name(), 3508, true));
  ils.push_back(boost::shared_ptr<Ils>(new IlsConvolution(d, ils_tab)));

  coeff = 4.74980283e+03, 1.99492886e-01;
  ils_tab.reset(new IlsTableLinear(hf, 2, "SC-Band", "strong_co2"));
  d.reset(new DispersionPolynomial(coeff, flag, units::inv_cm, 
				   ils_tab->band_name(), 2005, true));
  ils.push_back(boost::shared_ptr<Ils>(new IlsConvolution(d, ils_tab)));

  IlsInstrument inst(ils);
  StateVector sv;
  sv.add_observer(inst);
  Array<double,1> x(4);
  x = 1.28695614e+04, 5.74982835e+03,4.74980283e+03,0;
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
  SpectralDomain sd(wn_in, units::inv_cm);
  SpectralRange sr(rad_hres_in, units::dimensionless); // Ignore units
  BOOST_CHECK_MATRIX_CLOSE(inst.apply_instrument_model(Spectrum(sd, sr),
						       plist, 0).
			   spectral_range().data(), rad_out_expect);

  Array<double, 2> jac_rad_fake(rad_hres_in.rows(), 4);
  jac_rad_fake = 0;
  jac_rad_fake(Range::all(), 3) = rad_hres_in;
  ArrayAd<double, 1> rad_hres_in2(rad_hres_in, jac_rad_fake);
  SpectralRange sr2(rad_hres_in2, units::dimensionless); // Ignore units
  BOOST_CHECK_MATRIX_CLOSE(inst.apply_instrument_model(Spectrum(sd, sr), 
						       plist, 0).
			   spectral_range().data(), rad_out_expect);
  Array<double, 2> jac = 
    inst.apply_instrument_model(Spectrum(sd, sr2), plist, 0).
    spectral_range().data_ad().jacobian();
  Array<double, 1> v0 = 
    inst.apply_instrument_model(Spectrum(sd, sr2), plist, 0).
    spectral_range().data();
  double epsilon = 1e-3;
  x(0) += epsilon;
  sv.update_state(x);
  Array<double, 1> v1 = 
    inst.apply_instrument_model(Spectrum(sd, sr2), plist, 0).
    spectral_range().data();
  Array<double, 2> jacd(jac.shape());
  jacd = 0;
  jacd(Range::all(), 0) = (v1 - v0) / epsilon;
  x(0) -= epsilon;
  sv.update_state(x);
  rad_hres_in2.value() *= 1 + epsilon;
  SpectralRange sr3(rad_hres_in2, units::dimensionless); // Ignore units
  Array<double, 1> v2 = 
    inst.apply_instrument_model(Spectrum(sd, sr2), plist, 0).
    spectral_range().data();
  jacd(Range::all(), 3) = (v2 - v0) / epsilon;
  BOOST_CHECK_MATRIX_CLOSE(jac, jacd);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(ils_instrument_timing, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(timing)
{
  is_timing_test();
  const Instrument& inst = *config_instrument;
  IfstreamCs expected(test_data_dir() + "expected/ils_convolution/basic");
  Array<double, 1> wn_in, rad_hres_val;
  expected >> wn_in >> rad_hres_val;
  Array<double, 2> rad_hres_jac(rad_hres_val.rows(), 
				 config_state_vector->state().rows());
  for(int i = 0; i < rad_hres_jac.cols(); ++i)
    rad_hres_jac(Range::all(), i) = rad_hres_val;
  ArrayAd<double, 1> rad_hres(rad_hres_val, rad_hres_jac);
  std::vector<int> plist = 
    config_spectral_window->grid_indexes(inst.pixel_spectral_domain(0), 0);
  Spectrum s(SpectralDomain(wn_in, units::inv_cm),
	     SpectralRange(rad_hres, units::dimensionless));
  AccumulatedTimer tm("IlsInstrument");
  {
    FunctionTimer ft = tm.function_timer();
    Spectrum res = inst.apply_instrument_model(s, plist, 0);
  }
  std::cout << tm << "\n";
}

BOOST_AUTO_TEST_SUITE_END()

