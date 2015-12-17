begin FILE
  source_filename = oco_l2.inp

  ###############  
  # SOUNDING_INFO

  begin PICK
      section  = SOUNDING_INFO
      keyword  = spectrum_file
      template = <case_spectrum_filename>
  end PICK

end FILE
