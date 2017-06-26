# Modify noise column before adding noise to radiance
begin FILE
  source_filename = <case_spectrum_filename>

  begin MODIFY_DATA
     columns = 'noise'
     rows    = '0,1016'
     modify  = '<original> * <oco_noise_modify_1>'
     only_if = 'eval( len("<oco_noise_modify_1>") > 0 )'
     only_if = 'eval( "<instrument_name>".lower().find("oco") == 0 )'
  end MODIFY_DATA

  begin MODIFY_DATA
     columns = 'noise'
     rows    = '1016,2032'
     modify  = '<original> * <oco_noise_modify_2>'
     only_if = 'eval( len("<oco_noise_modify_2>") > 0 )'
     only_if = 'eval( "<instrument_name>".lower().find("oco") == 0 )'
  end MODIFY_DATA

  begin MODIFY_DATA
     columns = 'noise'
     rows    = '2032,3048'
     modify  = '<original> * <oco_noise_modify_3>'
     only_if = 'eval( len("<oco_noise_modify_3>") > 0 )'
     only_if = 'eval( "<instrument_name>".lower().find("oco") == 0 )'
  end MODIFY_DATA

  begin MODIFY_DATA
     columns = 'noise'
     modify  = '<original> * <gosat_noise_modify>'
     only_if = 'eval( len("<gosat_noise_modify>") > 0 )'
     only_if = 'eval( "<instrument_name>".lower() == "gosat" )'
  end MODIFY_DATA

end FILE

begin FILE
  source_filename = <case_spectrum_filename>

  # Add noise with cut off to bands seperately so
  # that mean noise value for each band is handled
  # seperately
  begin NOISIFY_SPECTRA_MODULE
     noise_cut_off = 2.0
     pixel_rows    = '0,1016'
     only_if       = 'eval( "<instrument_name>".lower().find("oco") == 0 )'
  end NOISIFY_SPECTRA_MODULE

  begin NOISIFY_SPECTRA_MODULE
     noise_cut_off = 2.0
     pixel_rows    = '1016,2032'
     only_if       = 'eval( "<instrument_name>".lower().find("oco") == 0 )'
  end NOISIFY_SPECTRA_MODULE

  begin NOISIFY_SPECTRA_MODULE
     noise_cut_off = 2.0
     pixel_rows    = '2032,3048'
     only_if       = 'eval( "<instrument_name>".lower().find("oco") == 0 )'
  end NOISIFY_SPECTRA_MODULE

  # Use a dumber approcach to adding noise for gosat
  # bands for now
  begin NOISIFY_SPECTRA_MODULE
     only_if       = 'eval( "<instrument_name>".lower() == "gosat" )'
  end NOISIFY_SPECTRA_MODULE

end FILE

