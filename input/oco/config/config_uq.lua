------------------------------------------------------------
--- Uncertainty Quantification Testing Configuration
------------------------------------------------------------

require "oco_base_config"

config = OcoBaseConfig:new()

require "uncertainty_quantification"
init_uq(config)

config:do_config()
