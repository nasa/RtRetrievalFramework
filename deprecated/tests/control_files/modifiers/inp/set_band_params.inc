# Set DISP retrievial indicies
begin FILE
  source_filename = oco_l2.inp

  begin PICK
      section  = PARAMETER_DEFINITION->INSTRUMENT
      keyword  = 'retrieval_indices(1)' 'retrieval_indices(2)' 'retrieval_indices(3)'
      template = ' '
      which_section = DISP
  end PICK

  begin PICK
      section  = PARAMETER_DEFINITION->BRDF->SPECTRALLY_DEPENDENT
      keyword  = 'retrieval_indices(1)' 'retrieval_indices(2)' 'retrieval_indices(3)'
      template = ' '
      only_if  = 'eval( "[brdf_type:case_sounding_id]" == "lambertian" )'
  end PICK


  begin MAP
     'ABO2' 1
     'WCO2' 2
     'SCO2' 3
  end MAP

  begin PICK
      section  = PARAMETER_DEFINITION->INSTRUMENT
      keyword  = 'retrieval_indices([band])'
      template = '1'
      which_section = DISP
      only_if  = 'eval( "<band>" == "ABO2" or "<band>" == "WCO2" or "<band>" == "SCO2" )'
  end PICK

  begin PICK
      section  = PARAMETER_DEFINITION->BRDF->SPECTRALLY_DEPENDENT
      keyword  = 'retrieval_indices([band])'
      template = '1 2'
      only_if  = 'eval( "[brdf_type:case_sounding_id]" == "lambertian" )' 
      only_if  = 'eval( "<band>" == "ABO2" or "<band>" == "WCO2" or "<band>" == "SCO2" )'
  end PICK

  begin PICK
      section  = PARAMETER_DEFINITION->INSTRUMENT
      keyword  = 'retrieval_indices(1)' 'retrieval_indices(2)' 
      template = '1'
      which_section = DISP
      only_if  = 'eval( "<band>" == "ABO2_WCO2" )'
  end PICK

  begin PICK
      section  = PARAMETER_DEFINITION->BRDF->SPECTRALLY_DEPENDENT
      keyword  = 'retrieval_indices(1)' 'retrieval_indices(2)'
      template = '1 2'
      only_if  = 'eval( "<band>" == "ABO2_WCO2" )'
  end PICK

  begin PICK
      section  = PARAMETER_DEFINITION->INSTRUMENT
      keyword  = 'retrieval_indices(1)' 'retrieval_indices(2)' 'retrieval_indices(3)'
      template = '1'
      which_section = DISP
      only_if  = 'eval( "<band>" == "3BND" )'
  end PICK

  begin PICK
      section  = PARAMETER_DEFINITION->BRDF->SPECTRALLY_DEPENDENT
      keyword  = 'retrieval_indices(1)' 'retrieval_indices(2)' 'retrieval_indices(3)'
      template = '1 2'
      only_if  = 'eval( "<band>" == "3BND" )'
  end PICK

end FILE
