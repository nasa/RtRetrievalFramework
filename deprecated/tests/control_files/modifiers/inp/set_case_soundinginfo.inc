begin FILE
  source_filename = oco_l2.inp

  ###############  
  # SOUNDING_INFO

  begin PICK
      section  = SOUNDING_INFO
      keyword  = soundinginfo_file
      template = <case_soundinginfo_filename>
  end PICK

end FILE