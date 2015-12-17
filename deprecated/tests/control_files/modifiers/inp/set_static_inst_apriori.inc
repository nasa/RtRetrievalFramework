begin FILE
  source_filename = oco_l2.inp

  ######
  # DISP

  begin MAP
     name = disp_ap_map
     begin VALUES
        oco        in/static/instrument/dispersion_orbit_sim.dat
        oco_tconn  in/static/instrument/dispersion_5.dat
        gosat      in/static/instrument/dispersion.dat
     end VALUES
  end MAP

  begin PICK
      section = PARAMETER_DEFINITION->INSTRUMENT
      keyword = a_priori
      template = [disp_ap_map:instrument_name]
      which_section = DISP
  end PICK

  #######
  # NOISE

  begin MAP
     name = noise_file_map
     begin VALUES
        oco       in/static/instrument/snr_coefs_5.dat
        oco_tconn in/static/instrument/snr_coefs_5.dat
        gosat     in/static/instrument/noise_cnv_P_H.dat
     end VALUES
  end MAP

  begin PICK
      section  = SOUNDING_INFO
      keyword  = noise_file
      template = [noise_file_map:instrument_name]
      only_if  = 'eval( "<use_run_mode>".upper() == "FORWARD_MODEL" or "<use_run_mode>".upper() == "JACOBIAN_ONLY" )'
  end PICK

  begin PICK
      section  = SOUNDING_INFO
      keyword  = noise_file
      delete   = True
      only_if  = 'eval( "<use_run_mode>".upper() != "FORWARD_MODEL" and "<use_run_mode>".upper() != "JACOBIAN_ONLY" )'
  end PICK

end FILE
