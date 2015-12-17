begin FILE
  source_filename = oco_l2.inp

  begin SETUP_ABSCO_FILES
      absco_path = <use_absco_path>
      prefer_file = 'with_lm'
      prefer_file = 'hi2008'
      prefer_file = '_alpha.abs'
      only_if     = 'eval( len("<use_absco_path>") > 0 )'
  end SETUP_ABSCO_FILES

end FILE
