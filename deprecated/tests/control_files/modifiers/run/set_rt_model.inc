begin FILE
  source_filename = oco_l2.run

  begin PICK
      section  = ALGORITHMS
      keyword  = rt_model_type
      template = <rt_type>
      only_if  = 'eval( len("<rt_type>") > 0 )'
  end PICK

end FILE
