------------------------------------------------------------
--- Default OCO configuration modified by Rail Road
--- Valley configuration.
------------------------------------------------------------

require "oco_base_config"

config = OcoBaseConfig:new()

require "rail_road_valley"
init_rrv(config)

-- Write jacobians and change convergence threshold
config.write_jacobian = true
config.solver.threshold = 1e-2

config:do_config()
