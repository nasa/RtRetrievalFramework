#include "atmosphere_fixture.h"
#include "pressure_fixed_level.h"
using namespace FullPhysics;
using namespace blitz;

AtmosphereFixture::AtmosphereFixture()
{
  atm = dynamic_cast<const AtmosphereOco&>(*config_atmosphere).clone();
  press_level = config_pressure_level_input;

  // Create a new SV for our cloned Atmosphere to use
  statev.reset(new StateVector);
  statev->add_observer(*atm);
  statev->update_state(config_initial_guess->initial_guess());
}

/// Set the surface pressure of the atmosphere.
void AtmosphereFixture::set_surface_pressure(double x)
{
  atm->set_surface_pressure_for_testing(x);
}

