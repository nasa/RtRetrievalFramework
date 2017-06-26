------------------------------------------------------------
--- This is a sample put together for Chris O'dell to change
--- the first guess based in an input initial guess file.
------------------------------------------------------------

require "oco_base_config"

config = OcoBaseConfig:new()

--- Change this to wherever the initial guess file is located.
local hfile_name = "/Users/smyth/Level2Build/build/initial_guess.h5"

-- Determine the column to use in the initial guess file. We get the
-- current directory, and strip off the ending _????, i.e, _0001. This is
-- a 1-based index, so we subtract by 1 to get the index.

local cdir = os.getenv("PWD")
-- This complicated expression just gets everything after the last '_'
local sindex = cdir:sub(cdir:len() - string.find(cdir:reverse(), '_') + 2)
local igindex = tonumber(sindex) - 1

local orig_fm_ig = config.fm.creator.initial_guess
config.fm.creator.initial_guess = function(self)
   local hfile = HdfFile(hfile_name)
   local nfg = hfile:read_double_2d("new_first_guess")
   local cfg = hfile:read_int_1d("change_first_guess")
   local ig = orig_fm_ig(self)
   local fg = ig:initial_guess()
   for i = 1,cfg:rows() do
      if(cfg(i-1) == 1) then
	 fg:set(i-1,nfg(i-1, igindex))
      end
   end
   ig = ConfigCommon.update_initial_guess(ig, fg)
   return ig
end

config:do_config()

