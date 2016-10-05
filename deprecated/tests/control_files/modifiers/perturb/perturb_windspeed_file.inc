begin FILE
   source_filename = <case_windspeed_filename>
   only_if  = 'eval( "<brdf_type>".lower() == "coxmunk" )'

   ###########
   # WINDSPEED
 
   # Perturb with random number
   begin MODIFY_DATA
      columns = WINDSPEED
      modify  = 'abs( <original> + random.gauss(0,5) )'
   end MODIFY_DATA
end FILE
