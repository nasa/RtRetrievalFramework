------------------------------------------------------------
--- Uses the orbit simulator base configuration instead
--- of the production base configuration.
------------------------------------------------------------

require "gosat_orbit_sim_base_config"

config = GosatOrbitSimBaseConfig:new()

config:do_config()
