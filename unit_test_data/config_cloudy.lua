------------------------------------------------------------
--- Do a run using the FixedLevelBaseConfig, with sounding changed to
--- one that is cloudy
------------------------------------------------------------

require "fixed_level_base_config"

config = FixedLevelBaseConfig:new {
   sid_string = "20090627205855"
}

config:do_config()
