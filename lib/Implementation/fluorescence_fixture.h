#include "fluorescence_effect.h"
#include "configuration_fixture.h"
#include <boost/shared_ptr.hpp>

using namespace FullPhysics;
using namespace blitz;

class FluorescenceFixture: public ConfigurationFixture {
public:
  FluorescenceFixture()
    : ConfigurationFixture("config_fluor.lua") 
  {
      epsilon.resizeAndPreserve(114);

      boost::shared_ptr<std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > > > vvse =
        lua_config["spectrum_effect"].value_ptr<std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > > >();
      config_fluor = boost::dynamic_pointer_cast<FluorescenceEffect>(vvse->at(0).at(1));
  }
  virtual ~FluorescenceFixture() {}
  boost::shared_ptr<FluorescenceEffect> config_fluor;
};
