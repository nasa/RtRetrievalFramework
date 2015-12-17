#ifndef AEROSOL_SHAPE_FIXTURE_H
#define AEROSOL_SHAPE_FIXTURE_H

#include "configuration_fixture.h"

namespace FullPhysics {

class AerosolShapeFixture: public ConfigurationFixture {
public:
  AerosolShapeFixture()
    : ConfigurationFixture("config_aer_shape.lua") {}
  virtual ~AerosolShapeFixture() {}
};

}

#endif
