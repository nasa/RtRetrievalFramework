begin KEYWORDS
   perturb_list = 'eval( "<perturb_list|T H2O CO2 O2>".upper().split() )'
end KEYWORDS

begin FILE
   source_filename = <case_atmosphere_filename>

   ##############
   # TEMPERATURE

   begin MODIFY_DATA
      columns     = T
      modify      = 'abs( <original> + <temp_offset> )'

      begin CONSTANTS
         temp_offset = 'eval( random.gauss(0,5))'
      end CONSTANTS

      only_if   = 'eval( "T" in "<perturb_list>" )'
   end MODIFY_DATA

   #####
   # H2O
   begin MODIFY_DATA
      columns   = H2O
      modify    = 'abs( <original> * ( 1 + <h2o_scale> ) )'

      begin CONSTANTS
         h2o_scale = 'eval( random.gauss(0,0.316) )'
      end CONSTANTS

      only_if   = 'eval( "H2O" in "<perturb_list>" )'
   end MODIFY_DATA

   #####
   # CO2
   begin COVARIANCE_REALIZATION
      column     = CO2
      covariance = in/static/covariance/co2_cov_20.dat

      only_if   = 'eval( "CO2" in "<perturb_list>" )'
   end COVARIANCE_REALIZATION

   #####
   # O2
   begin MODIFY_DATA
      columns   = O2
      modify    = '.2095'
      add_column = True

      only_if   = 'eval( "O2" in "<perturb_list>" )'
   end MODIFY_DATA

end FILE

begin FILE
   source_filename = <case_atmosphere_filename>

   begin MODIFY_DATA 
      columns = T
      delete  = true 
      not_if   = 'eval( "T" in "<perturb_list>" )' 
   end MODIFY_DATA 

   begin MODIFY_DATA
      columns = H2O
      delete  = true
      not_if   = 'eval( "H2O" in "<perturb_list>" )'
   end MODIFY_DATA

   begin MODIFY_DATA
      columns = CO2
      delete  = true
      not_if   = 'eval( "CO2" in "<perturb_list>" )'
   end MODIFY_DATA

   begin MODIFY_DATA
      columns = O2
      delete  = true
      not_if   = 'eval( "O2" in "<perturb_list>" )'
   end MODIFY_DATA

end FILE
