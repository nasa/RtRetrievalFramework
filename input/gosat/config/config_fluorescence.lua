------------------------------------------------------------
--- Adds the use of FluorescenceEffect to the retrieval
------------------------------------------------------------

require "gosat_base_config"

config = GosatBaseConfig:new()

table.insert(config.fm.spectrum_effect.speceff, "fluorescence")

config:do_config()
