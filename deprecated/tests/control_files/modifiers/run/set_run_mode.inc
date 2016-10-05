begin FILE
  source_filename = oco_l2.run

  begin PICK
      section  = CONTROL
      keyword  = run_mode
      template = <use_run_mode>
  end PICK

  begin PICK
      section  = OUTPUT
      keyword  = high_res_spectra
      template = TRUE
      only_if  = 'eval( "<run_mode>".upper() == "FORWARD_MODEL" )'
  end PICK

end FILE
