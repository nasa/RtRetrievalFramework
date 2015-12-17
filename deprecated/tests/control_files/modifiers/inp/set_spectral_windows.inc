begin FILE
  source_filename = oco_l2.inp

  begin PICK
      section  = WINDOW_INFO
      keyword  = spectral_windows(1)
      template = '<spectral_windows_1>'
      only_if  = 'eval(len("<spectral_windows_1>") > 0)'
  end PICK

  begin PICK
      section  = WINDOW_INFO
      keyword  = spectral_windows(2)
      template = '<spectral_windows_2>'
      only_if  = 'eval(len("<spectral_windows_2>") > 0)'
  end PICK

  begin PICK
      section  = WINDOW_INFO
      keyword  = spectral_windows(3)
      template = '<spectral_windows_3>'
      only_if  = 'eval(len("<spectral_windows_3>") > 0)'
  end PICK

end FILE