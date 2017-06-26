#include "unit_test_support.h"
#include "aerosol_param_output.h"
#include "aerosol_shape_fixture.h"
#include "output_hdf.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(aerosol_param_output, AerosolShapeFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<AerosolOptical> a =
    boost::dynamic_pointer_cast<AerosolOptical>(config_aerosol);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("aerosol_param_output.h5", config_pressure->number_level(), config_state_vector->state().rows(), config_aerosol->number_particle(), 3));
  add_file_to_cleanup("aerosol_param_output.h5");
  for(int aer_idx = 0; aer_idx < config_aerosol->number_particle(); aer_idx++) {
    boost::shared_ptr<AerosolExtinctionImpBase> aer_ext = 
      boost::dynamic_pointer_cast<AerosolExtinctionImpBase>(a->aerosol_extinction(aer_idx));
    AerosolParamOutput po(aer_ext);
    po.register_output_apriori(out);
    po.register_output(out);
  }

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in aerosol unit test.
  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


