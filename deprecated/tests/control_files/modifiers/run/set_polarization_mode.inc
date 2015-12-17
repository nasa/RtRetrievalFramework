begin FILE
  source_filename = oco_l2.run

  begin PICK
      section  = ALGORITHMS
      keyword  = polarization
      template = '<use_polarization>'
      only_if  = 'eval( len("<use_polarization>") > 0 )'
  end PICK

end FILE
