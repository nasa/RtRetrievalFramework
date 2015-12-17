begin FILE
  source_filename = oco_l2.run

  begin PICK
      section  = ALGORITHMS
      keyword  = solar_doppler_shift
      template = <use_solar_doppler>
      only_if  = 'eval( len("<use_solar_doppler>") > 0 )'
  end PICK

  begin PICK
      section  = ALGORITHMS
      keyword  = inst_doppler_shift
      template = <use_inst_doppler>
      only_if  = 'eval( len("<use_inst_doppler>") > 0 )'
  end PICK

end FILE
