begin FILE
  source_filename = oco_l2.inp

  #########
  # AEROSOL

  begin PICK
      section = PARAMETER_DEFINITION->AEROSOL
      keyword = a_priori
      template = <case_aerosol_filename>
  end PICK

end FILE


