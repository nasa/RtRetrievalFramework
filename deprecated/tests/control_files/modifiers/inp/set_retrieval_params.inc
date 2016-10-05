# Removes any retrieval indicies, scale parameters not in the desired retrieval vector

begin FILE
  source_filename = oco_l2.inp

  ##################
  # retrieval_vector

  # Set new retrieval vector specified by user
  begin PICK
     section  = PARAMETER_DEFINITION
     keyword  = retrieval_vector
     template = <use_retrieval_vector>
     only_if  = 'eval( len("<use_retrieval_vector>") > 0 )'
  end PICK

  #####
  # CO2

  begin PICK
      section  = PARAMETER_DEFINITION->GAS
      keyword  = retrieval_indices scale_parameters
      delete   = True
      which_section = CO2
      not_if   = 'eval( "CO2" in "<use_retrieval_vector>".upper().split() )'
      only_if  = 'eval( len("<use_retrieval_vector>") > 0 )'
  end PICK

  #####
  # H2O

  begin PICK
      section  = PARAMETER_DEFINITION->GAS
      keyword  = retrieval_indices scale_parameters
      delete   = True
      which_section = H2O
      not_if   = 'eval( "H2O" in "<use_retrieval_vector>".upper().split() )'
      only_if  = 'eval( len("<use_retrieval_vector>") > 0 )'
  end PICK

  ####
  # O2

  begin PICK
      section  = PARAMETER_DEFINITION->GAS
      keyword  = retrieval_indices scale_parameters
      delete   = True
      which_section = O2
      not_if   = 'eval( "O2" in "<use_retrieval_vector>".upper().split() )'
      only_if  = 'eval( len("<use_retrieval_vector>") > 0 )'
  end PICK

  #######
  # PSURF

  begin PICK
      section  = PARAMETER_DEFINITION->SURFACE_PRESSURE
      keyword  = retrieval_indices 
      delete   = True
      not_if   = 'eval( "PSURF" in "<use_retrieval_vector>".upper().split() )'
      only_if  = 'eval( len("<use_retrieval_vector>") > 0 )'
  end PICK

  ###
  # T

  begin PICK
      section  = PARAMETER_DEFINITION->TEMPERATURE
      keyword  = retrieval_indices scale_parameters
      delete   = True
      not_if   = 'eval( "T" in "<use_retrieval_vector>".upper().split() )'
      only_if  = 'eval( len("<use_retrieval_vector>") > 0 )'
  end PICK

  ########
  # ALBEDO

  begin PICK
      section  = PARAMETER_DEFINITION->BRDF->SPECTRALLY_DEPENDENT
      keyword  = retrieval_indices(1) retrieval_indices(2) retrieval_indices(3)
      delete   = True
      not_if   = 'eval( "LAMBERTIAN" in "<use_retrieval_vector>".upper().split() )'
      only_if  = 'eval( "<brdf_type>".lower() == "lambertian" )'
      only_if  = 'eval( len("<use_retrieval_vector>") > 0 )'
  end PICK

  ###########
  # WINDSPEED

  begin PICK
      section  = PARAMETER_DEFINITION->BRDF->SPECTRALLY_INDEPENDENT
      keyword  = retrieval_indices(1) retrieval_indices(2) retrieval_indices(3)
      delete   = True
      not_if   = 'eval( "COXMUNK" in "<use_retrieval_vector>".upper().split() )'
      only_if  = 'eval( "<brdf_type>".lower() == "coxmunk" )'
      only_if  = 'eval( len("<use_retrieval_vector>") > 0 )'
  end PICK

  ######
  # DISP

  begin PICK
      section  = PARAMETER_DEFINITION->INSTRUMENT
      keyword  = retrieval_indices(1) retrieval_indices(2) retrieval_indices(3)
      delete   = True
      not_if   = 'eval( "DISP" in "<use_retrieval_vector>".upper().split() )'
      only_if  = 'eval( len("<use_retrieval_vector>") > 0 )'
      which_section = DISP
  end PICK

  #########
  # AEROSOL

  begin PICK
      section  = PARAMETER_DEFINITION->AEROSOL
      keyword  = retrieval_indices(1) retrieval_indices(2) retrieval_indices(3)
      delete   = True
      ignore_missing = True # Ignore missing 
      which_section = 'eval( "<use_retrieval_vector>".upper().split() )'
      only_if  = 'eval( len("<use_retrieval_vector>") > 0 )'
  end PICK

end FILE
