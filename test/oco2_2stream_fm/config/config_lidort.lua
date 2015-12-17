-- No changes to base needed to use LIDORT as it is the default

dofile "base_test.lua"

config.fm.rt.lsi_constant.dedicated_twostream = false 

config:do_config()
