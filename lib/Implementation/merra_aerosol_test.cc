#include "merra_aerosol.h"
#include "aerosol_optical.h"
#include "pressure_sigma.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(merra_aerosol, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
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
  boost::shared_ptr<HdfFile> mclimate;
  mclimate.reset(new HdfFile(merra_data_dir() + "MERRA_Composite_Selection_M02_O2A.hdf5"));
  Array<double, 2> cov(3, 3);
  cov =
    3.24, 0, 0,
    0, 0.16, 0,
    0, 0, 1e-4;
  MerraAerosol ma(*mclimate, aprop, DoubleWithUnit(36.278789520264, "deg"),
		  DoubleWithUnit(140.2403717041, "deg"), p, rh_dummy, cov);
  boost::shared_ptr<AerosolOptical> aop = boost::dynamic_pointer_cast<AerosolOptical>(ma.aerosol());
  BOOST_CHECK(aop->number_particle() == 2);
  BOOST_CHECK(aop->aerosol_name()[0] == "DU");
  BOOST_CHECK(aop->aerosol_name()[1] == "SO");
  boost::shared_ptr<CompositeInitialGuess> ig =
    boost::dynamic_pointer_cast<CompositeInitialGuess>(ma.initial_guess());
  Array<double, 1> ig_expect(6);
  ig_expect =  -2.95423, 0.548065, 0.185876, -3.06037, 0.92695, 0.378612;
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

BOOST_AUTO_TEST_CASE(dateline)
{
  // Check correct handling near the dateline

  blitz::Array<double, 1> a(3), b(3);
  a = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  boost::shared_ptr<Pressure> p (new PressureSigma(a,b, psurf, true));
  boost::shared_ptr<RelativeHumidity> rh_dummy;
  HdfFile aprop(test_data_dir() + "l2_merra_aerosol.h5");
  boost::shared_ptr<HdfFile> mclimate;
  mclimate.reset(new HdfFile(merra_data_dir() + "MERRA_Composite_Selection_M02_O2A.hdf5"));
  Array<double, 2> cov(3, 3);
  cov =
    3.24, 0, 0,
    0, 0.16, 0,
    0, 0, 1e-4;
  MerraAerosol ma(*mclimate, aprop, DoubleWithUnit(36.278789520264, "deg"),
		  DoubleWithUnit(179.4, "deg"), p, rh_dummy, cov);
  boost::shared_ptr<AerosolOptical> aop = boost::dynamic_pointer_cast<AerosolOptical>(ma.aerosol());
  BOOST_CHECK(aop->number_particle() == 2);
  BOOST_CHECK(aop->aerosol_name()[0] == "SS");
  BOOST_CHECK(aop->aerosol_name()[1] == "DU");
}

BOOST_AUTO_TEST_CASE(dateline2)
{
  // Check correct handling near the dateline

  blitz::Array<double, 1> a(3), b(3);
  a = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  boost::shared_ptr<Pressure> p (new PressureSigma(a,b, psurf, true));
  boost::shared_ptr<RelativeHumidity> rh_dummy;
  HdfFile aprop(test_data_dir() + "l2_merra_aerosol.h5");
  boost::shared_ptr<HdfFile> mclimate;
  mclimate.reset(new HdfFile(merra_data_dir() + "MERRA_Composite_Selection_M02_O2A.hdf5"));
  Array<double, 2> cov(3, 3);
  cov =
    3.24, 0, 0,
    0, 0.16, 0,
    0, 0, 1e-4;
  MerraAerosol ma(*mclimate, aprop, DoubleWithUnit(36.278789520264, "deg"),
		  DoubleWithUnit(-180.001, "deg"), p, rh_dummy, cov);
  boost::shared_ptr<AerosolOptical> aop = boost::dynamic_pointer_cast<AerosolOptical>(ma.aerosol());
  BOOST_CHECK(aop->number_particle() == 2);
  BOOST_CHECK(aop->aerosol_name()[0] == "SS");
  BOOST_CHECK(aop->aerosol_name()[1] == "DU");
}

BOOST_AUTO_TEST_CASE(dejian_example)
{
  // These are some examples coming directly from Dejian, this makes
  // sure we are doing the calculation the same way.

  // Dummy pressure. We need to pass this in to create the Aerosol,
  // but we don't actually use this in this test so any pressure is ok.
  blitz::Array<double, 1> a(3), b(3);
  a = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  boost::shared_ptr<Pressure> p (new PressureSigma(a,b, psurf, true));
  boost::shared_ptr<RelativeHumidity> rh_dummy;
  HdfFile aprop(test_data_dir() + "l2_merra_aerosol.h5");
  boost::shared_ptr<HdfFile> mclimate;
  mclimate.reset(new HdfFile(merra_data_dir() + "MERRA_Composite_Selection_M01_O2A.hdf5"));
  Array<double, 2> cov(3, 3);
  cov =
    3.24, 0, 0,
    0, 0.16, 0,
    0, 0, 1e-4;
  MerraAerosol ma1(*mclimate, aprop, DoubleWithUnit(17.000, "deg"),
		   DoubleWithUnit(8, "deg"), p, rh_dummy, cov);
  MerraAerosol ma2(*mclimate, aprop, DoubleWithUnit(36.6040, "deg"),
		   DoubleWithUnit(-97.4860, "deg"), p, rh_dummy, cov);
  MerraAerosol ma3(*mclimate, aprop, DoubleWithUnit(39.9060, "deg"),
		   DoubleWithUnit(116.4060, "deg"), p, rh_dummy, cov);
  MerraAerosol ma4(*mclimate, aprop, DoubleWithUnit(-7.9165, "deg"),
		   DoubleWithUnit(-14.3325, "deg"), p, rh_dummy, cov);
  boost::shared_ptr<AerosolOptical> aop1 = boost::dynamic_pointer_cast<AerosolOptical>(ma1.aerosol());
  boost::shared_ptr<AerosolOptical> aop2 = boost::dynamic_pointer_cast<AerosolOptical>(ma2.aerosol());
  boost::shared_ptr<AerosolOptical> aop3 = boost::dynamic_pointer_cast<AerosolOptical>(ma3.aerosol());
  boost::shared_ptr<AerosolOptical> aop4 = boost::dynamic_pointer_cast<AerosolOptical>(ma4.aerosol());
  BOOST_CHECK(aop1->number_particle() == 2);
  BOOST_CHECK(aop1->aerosol_name()[0] == "DU");
  BOOST_CHECK(aop1->aerosol_name()[1] == "SO");
  BOOST_CHECK(aop2->number_particle() == 2);
  BOOST_CHECK(aop2->aerosol_name()[0] == "SO");
  BOOST_CHECK(aop2->aerosol_name()[1] == "OC");
  BOOST_CHECK(aop3->number_particle() == 2);
  BOOST_CHECK(aop3->aerosol_name()[0] == "DU");
  BOOST_CHECK(aop3->aerosol_name()[1] == "SO");
  BOOST_CHECK(aop4->number_particle() == 2);
  BOOST_CHECK(aop4->aerosol_name()[0] == "SS");
  BOOST_CHECK(aop4->aerosol_name()[1] == "OC");
}

BOOST_AUTO_TEST_SUITE_END()


