begin FILE
  source_filename = oco_l2.run

  begin PICK
      section  = ALGORITHMS
      keyword  = streams
      template = 8
      not_if   = '<use_lsi|True>'
  end PICK

  begin PICK
      section  = ALGORITHMS
      keyword  = lsi_correction
      template = False
      not_if   = '<use_lsi|True>'
  end PICK

end FILE
