begin SETTINGS
   # Define pressure levels file to use for each run, will be used in resampling
   pressure_levels_filename = in/static/scene/atmosphere/pressure_levels_20.dat

   # Used when these are missing
   aerosol_perturb_filename = in/scene/aerosol/aerosol_perturb.dat
   aerosol_covariance_filename = in/scene/aerosol/aerosol_cov.dat

   # To match the spectra created by the orbit simulator
   use_solar_linelist = 'solar_di_20100208.101'
   use_solar_cont = 'solar_cont_v2.dat'
   use_absco_path = '/groups/algorithm/l2_fp/absco/v3.1.0_alpha/'
   use_solar_doppler = False
   use_inst_doppler = False
   disp_index_scheme = 0

end SETTINGS
