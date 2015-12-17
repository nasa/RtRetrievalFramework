#include "ils_fts.h"
#include "level_1b_fts.h"
#include "spectral_window_range.h"
#include "state_vector.h"
#include "unit_test_support.h"
#include "radiance_scaling_sv_fit.h"
#include "dispersion_polynomial.h"
#include "ils_instrument.h"
#include <iostream>
using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ils_fts, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::vector<std::string> spectra;
  //one spectrometers
  spectra.push_back(test_data_dir() + 
		    "in/l1b/spec/pa20091103saaaaa_100223160344.008");

  ArrayWithUnit<double, 2> spec_range;
  spec_range.value.resize(1, 2);
  spec_range.units = units::inv_cm;
  spec_range.value = 6173.0, 6275.0;

  boost::shared_ptr<Level1bFts> l1b_fts(new Level1bFts(test_data_dir() + "in/pa20091103_100223163011.grl", spectra, spec_range));

  // Set up dispersion values
  blitz::Array<double, 1> disp_coeff(2);

  // Use value calculated (wrongly) by old code, the offset is not correct
  // and certainly not disvisible by the spacing which is correct
  disp_coeff = 6173.0156278633530746, 0.0075330826756407334374; 
  
  Array<bool, 1> flag_disp(2);
  Array<double, 2> disp_perturb(1, 2);
  flag_disp = true, false;
  disp_perturb = 0.01;
  boost::shared_ptr<DispersionPolynomial> 
    disp(new DispersionPolynomial(disp_coeff, flag_disp, "cm^-1", "Band 1", 
				   l1b_fts->radiance(0).data().rows(), false));
  
  boost::shared_ptr<Ils> ils_fts(new IlsFts(disp, disp_perturb, l1b_fts, 0, 
					    "Band 1", "band_1"));

  //// Test case we captured
  IfstreamCs indata(test_data_dir() + "/expected/ils_fts/ils_fts");
  Array<double, 1> wn_in, rad_hres_in, ils_out_expect, rad_out_expect;
  indata >> wn_in >> rad_hres_in >> ils_out_expect >> rad_out_expect;

  std::vector<int> pixel_list;
  for(int i = 10; i <= 100; ++i)
    pixel_list.push_back(i);
  SpectralDomain sd(wn_in, units::inv_cm);
  SpectralRange sr(rad_hres_in, units::dimensionless); // Ignore units
  Array<double,1> ils_conv = 
    ils_fts->apply_ils(wn_in, rad_hres_in, pixel_list);

  BOOST_CHECK_MATRIX_CLOSE_TOL(ils_out_expect, ils_conv, 2e-5);

  // Now create a full instrument, and compare with the full radiance
  // output

  // Set up continuum values
  Array<bool, 1> flag_cont(2);
  flag_cont = false;

  blitz::Array<double, 1> continuum(2);
  continuum = 2.03662741E+00, -1.16125220E-03;
  boost::shared_ptr<InstrumentCorrection> 
    ic(new RadianceScalingSvFit(continuum, flag_cont, 
			   DoubleWithUnit(0.5 * (6173.0 + 6275.0), "cm^-1"),
			   "Band 1"));
  std::vector<boost::shared_ptr<InstrumentCorrection> > ic2;
  ic2.push_back(ic);
  std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> >  > iclist;
  iclist.push_back(ic2);
  std::vector<boost::shared_ptr<Ils> > ils;
  ils.push_back(ils_fts);

  IlsInstrument inst_fts(ils, iclist);
  Array<double,1> rad_conv = 
    inst_fts.apply_instrument_model(Spectrum(sd, sr), pixel_list, 0).
    spectral_range().data();

  BOOST_CHECK_MATRIX_CLOSE_TOL(rad_out_expect(Range(10, 100)), rad_conv, 2e-5);
  
}

BOOST_AUTO_TEST_SUITE_END()
