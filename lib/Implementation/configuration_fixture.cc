#include "configuration_fixture.h"
#include <fenv.h>

using namespace FullPhysics;
using namespace blitz;

std::map<std::string, boost::shared_ptr<LuaState> > ConfigurationFixture::config;
ConfigurationFixture::ConfigurationFixture(const std::string& Config_file)
  : epsilon(112)
{
  if(!config[Config_file]) {
    // Disable floating point exeptions while loading 
    // Lua configuration due to the way Lua parses certain
    // things.

    // Mac doesn't have this function, even though it is a C99
    // function. We check for this during configuration.
#ifdef HAVE_FEENABLEEXCEPT
    fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
    config[Config_file] = LuaState::load_file(test_data_dir() + Config_file);
#ifdef HAVE_FEENABLEEXCEPT
    feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
  }

  lua_config = config[Config_file]->globals()["config"];
  epsilon = 1e-6;		  // Default
  epsilon(Range(0, 19)) = 1e-7;	  // CO2 VMR
  epsilon(21) = 1e-3;		  // Surface Pressure
  epsilon(22) = 1e-4;		  // Temperature
  epsilon(Range(23, 102)) = 1e-6; // Aerosol optical depth

  config_absorber = lua_config["absorber"].value_ptr<Absorber>();
  // Allow this to fail, we don't have aerosols if we happen to have a 
  // Rayleigh only atmosphere
  try {
    config_aerosol = lua_config["aerosol"].value_ptr<Aerosol>();
  } catch(const std::exception& e) {
    ;
  }
  config_atmosphere = lua_config["atmosphere"].value_ptr<RtAtmosphere>();
  config_state_vector = lua_config["state_vector"].value_ptr<StateVector>();
  if (!lua_config["pinp"].is_nil()) // Not all configs use this
    config_pressure_level_input = 
      lua_config["pinp"].value_ptr<PressureLevelInput>();
  config_pressure = lua_config["pressure"].value_ptr<Pressure>();
  config_register_output = 
    lua_config["register_output"].
    value<std::vector<boost::shared_ptr<RegisterOutputBase> > >();
  config_instrument = lua_config["instrument"].value_ptr<Instrument>();
  config_spectral_window = lua_config["spec_win"].value_ptr<SpectralWindow>();
  config_initial_guess = 
    lua_config["initial_guess"].value_ptr<InitialGuess>();
  config_solver = lua_config["conn_solver"].value_ptr<ConnorSolver>();
  // Allow this to fail, we don't have ground if we happen to be looking
  // up.
  try {
    config_ground = lua_config["ground"].value_ptr<Ground>();
  } catch(const std::exception& e) {
    ;
  }
  config_temperature = lua_config["temperature"].value_ptr<Temperature>();
  config_spectrum_sampling = 
    lua_config["spec_samp"].value_ptr<SpectrumSampling>();
  config_error_analysis = 
    lua_config["error_analysis"].value_ptr<ErrorAnalysis>();
  config_level_1b = lua_config["l1b"].value_ptr<Level1b>();
  config_rt = lua_config["rt"].value_ptr<RadiativeTransfer>();
  config_forward_model = lua_config["forward_model"].value_ptr<ForwardModel>();
  sv_initial.reference(config_initial_guess->initial_guess());
  config_state_vector->update_state(sv_initial);
}

ConfigurationCoxmunkFixture::ConfigurationCoxmunkFixture(const std::string& Config_file) 
  : ConfigurationFixture(Config_file)
{
  // Shrink array from base constructor less elements but first many should be the same
  epsilon.resizeAndPreserve(110);
}

