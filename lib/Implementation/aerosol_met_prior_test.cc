#include "aerosol_met_prior.h"
#include "aerosol_optical.h"
#include "pressure_sigma.h"
#include "unit_test_support.h"
#include "oco_sounding_id.h"
#include <boost/make_shared.hpp>

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(aerosol_met_prior, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  return; 			// Depends on data we don't normally
				// have available
  // Dummy pressure. We need to pass this in to create the Aerosol,
  // but we don't actually use this in this test so any pressure is ok.
  blitz::Array<double, 1> a(3), b(3);
  a = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  boost::shared_ptr<Pressure> p (new PressureSigma(a,b, psurf, true));
  // Dummy RelativeHumidity, we need to create Aerosol but not
  // actually used.
  boost::shared_ptr<RelativeHumidity> rh_dummy;
  HdfFile aprop(test_data_dir() + "l2_merra_aerosol.h5");
  // Temp file name
  std::string fname = "/home/smyth/Local/Level2/oco2_L2MetGL_05297a_150701_Bxxxx_170127155155d.h5";
  HdfFile met_hfile(fname);
  boost::shared_ptr<HdfSoundingId> sid =
    boost::make_shared<OcoSoundingId>(met_hfile, "2015070100334532");
  OcoMetFile met_file(fname, sid);
  Array<double, 2> cov(3, 3);
  cov =
    3.24, 0, 0,
    0, 0.16, 0,
    0, 0, 1e-4;
  AerosolMetPrior ma(met_file, aprop, p, rh_dummy, cov);
  boost::shared_ptr<AerosolOptical> aop = boost::dynamic_pointer_cast<AerosolOptical>(ma.aerosol());
  std::cerr << *aop << "\n";
  BOOST_CHECK(aop->number_particle() == 2);
  BOOST_CHECK(aop->aerosol_name()[0] == "SS");
  BOOST_CHECK(aop->aerosol_name()[1] == "SO");
  boost::shared_ptr<CompositeInitialGuess> ig =
    boost::dynamic_pointer_cast<CompositeInitialGuess>(ma.initial_guess());
  Array<double, 1> ig_expect(6);
  ig_expect =  -2.18912, 0.927757, 0.0723657, -5.54738, 0.691174, 0.261427;
  Array<double, 2> cov_expect(6,6);
  cov_expect =     
    3.24, 0, 0, 0, 0, 0,
    0, 0.16, 0, 0, 0, 0,
    0, 0, 1e-4, 0, 0, 0,
    0, 0, 0, 3.24, 0, 0,
    0, 0, 0, 0, 0.16, 0,
    0, 0, 0, 0, 0, 1e-4;
  BOOST_CHECK_MATRIX_CLOSE_TOL(ig->initial_guess(), ig_expect, 1e-4);
  BOOST_CHECK_MATRIX_CLOSE(ig->apriori_covariance(), cov_expect);
}

BOOST_AUTO_TEST_SUITE_END()


