begin KEYWORDS
   aerosol_types_used = '<aerosol_types|AERO1 AERO2 WC IC>'

   aerosol_apriori_type = '<aerosol_apriori_type|profile>'
   aerosol_shape_types  = '<aerosol_shape_types|gauss_4_types_lin>'
end KEYWORDS

# Delete file from scenarios since not using it, deleted
# to remove confusion
begin CLEANUP
   fail_on_error = False

   begin DELETE
      <case_aerosol_filename>
   end DELETE
end CLEANUP

begin FILE
  source_filename = <case_aod_filename>

  begin RANGE_MAP
     identifier = 'aod_scaling'

     begin RANGES
        '005' '0.0, 0.1'
        '015' '0.1, 0.2'
        '025' '0.2, 9999.0'
     end RANGES
  end RANGE_MAP
end FILE

begin FILE
  source_filename = oco_l2.inp

  begin MAP
      name = aer_shape_cov
      begin VALUES
         gauss_4_types_lin './in/static/covariance/aerosol_gauss_lin_cov.dat'
         gauss_4_types_log './in/static/covariance/aerosol_gauss_log_cov.dat'
         gauss_3_types_lin './in/static/covariance/aerosol_gauss_lin_cov.dat'
         gauss_3_types_log './in/static/covariance/aerosol_gauss_log_cov.dat'
      end VALUES
  end MAP

  # Remove all blocks and types if scaling of 0.00
  begin PICK
      section = PARAMETER_DEFINITION->AEROSOL
      delete  = True
      only_if = 'eval( "<aod_scaling>" == "000" )'
  end PICK

  begin PICK
      section  = PARAMETER_DEFINITION
      keyword  = aerosol_types
      template = ''
      only_if  = 'eval( "<aod_scaling>" == "000" )'
  end PICK

  # Set correct apriori profile for aerosol
  begin PICK
      section  = PARAMETER_DEFINITION->AEROSOL
      keyword  = a_priori
      template =  ./in/static/scene/aerosol/profiles/aerosol_<aod_scaling>_log_20.dat
      not_if   = 'eval( "<aod_scaling>" == "000" )'
      only_if  = 'eval( "<aerosol_apriori_type>".lower() == "profile" )' 
  end PICK

  # Set correct apriori profile for aerosol
  begin PICK
      section  = PARAMETER_DEFINITION->AEROSOL
      keyword  = a_priori
      template =  ./in/static/scene/aerosol/shapes/aerosol_<aerosol_shape_types>_<aod_scaling>.dat
      not_if   = 'eval( "<aod_scaling>" == "000" )'
      only_if  = 'eval( "<aerosol_apriori_type>".lower() == "shape" )' 
  end PICK

  begin PICK
      section  = PARAMETER_DEFINITION->AEROSOL
      keyword  = covariance
      template = [aer_shape_cov:aerosol_shape_types]
      not_if   = 'eval( "<aod_scaling>" == "000" )'
      only_if  = 'eval( "<aerosol_apriori_type>".lower() == "shape" )' 
  end PICK

  begin CREATE_AEROSOL_BLOCKS
      profile_names          = AERO1
      type_names             = vij_5
      set_retrieval_vector   = true
      not_if   = 'eval( "<aod_scaling>" == "000" )'
      only_if  = 'eval( "AERO1" in "<aerosol_types_used>".split() >= 0 )'
  end CREATE_AEROSOL_BLOCKS

  begin CREATE_AEROSOL_BLOCKS
      profile_names          = AERO2
      type_names             = vij_7
      set_retrieval_vector   = true
      not_if   = 'eval( "<aod_scaling>" == "000" )'
      only_if  = 'eval( "AERO2" in "<aerosol_types_used>".split() >= 0 )'
  end CREATE_AEROSOL_BLOCKS

  begin CREATE_AEROSOL_BLOCKS
      profile_names          = WC
      type_names             = wc_008
      set_retrieval_vector   = true
      not_if   = 'eval( "<aod_scaling>" == "000" )'
      only_if  = 'eval( "WC" in "<aerosol_types_used>".split() >= 0 )'
      not_if   = 'eval( "<aerosol_shape_types>".find("gauss_3_types") >= 0 )'
  end CREATE_AEROSOL_BLOCKS

  begin CREATE_AEROSOL_BLOCKS
      profile_names          = IC
      type_names             = ic_070
      set_retrieval_vector   = true
      not_if   = 'eval( "<aod_scaling>" == "000" )'
      only_if  = 'eval( "IC" in "<aerosol_types_used>".split() >= 0 )'
  end CREATE_AEROSOL_BLOCKS

end FILE
