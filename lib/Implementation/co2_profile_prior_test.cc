#include "co2_profile_prior.h"
#include "unit_test_support.h"
#include "oco_sounding_id.h"
#include "pressure_sigma.h"
#include <boost/make_shared.hpp>

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(co2_profile_prior, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  // return; 			// Depends on data we don't normally
				// have available
  // Temp file name
  std::string met_fname = "/home/smyth/Local/Level2/oco2_L2MetGL_05273a_150629_B9206p01_190718225702s.h5";
  std::string pr_fname = "/home/smyth/Local/Level2/oco2_L2CPrGL_05273a_150629_B9206p02_190719204712s.h5";
  HdfFile met_hfile(met_fname);
  boost::shared_ptr<HdfSoundingId> sid =
    boost::make_shared<OcoSoundingId>(met_hfile, "2015062909003032");
  OcoMetFile met(met_fname, sid);
  
  HdfFile hdf_static_input(test_data_dir() + "../input/oco/input/l2_oco_static_input.h5");
  blitz::Array<double, 1> sigma_a = hdf_static_input.read_field<double, 1>("Pressure/Pressure_sigma_a");
  blitz::Array<double, 1> sigma_b = hdf_static_input.read_field<double, 1>("Pressure/Pressure_sigma_b");
  double surface_pressure = met.surface_pressure();
  boost::shared_ptr<PressureSigma> press
    (new PressureSigma(sigma_a, sigma_b, surface_pressure, false));
  std::cerr << *press;
}

BOOST_AUTO_TEST_SUITE_END()


