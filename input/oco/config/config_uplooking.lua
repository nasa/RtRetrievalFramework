------------------------------------------------------------
-- Makes the  necessary changes to the nominal OCO
-- configuration to support uplooking data retrievals.
------------------------------------------------------------

require "oco_uplooking_tvac_base_config"

tvac_config = OcoUplookingBaseConfig:new()

tvac_config:do_config()
