-- Add option to enable 2stream as RT used in LSI 

dofile "base_test.lua"

config.fm.rt.lsi_constant.dedicated_twostream = true
config.fm.rt.lsi_constant.twostream_full_quadrature = false

config:do_config()
