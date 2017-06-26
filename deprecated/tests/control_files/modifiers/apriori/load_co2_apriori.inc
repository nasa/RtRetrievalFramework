begin KEYWORDS
   co2_apriori_filename   = '<co2_apriori_filename|in/scene/atmosphere/co2_apriori.dat>'
   num_co2_apriori_levels = '<num_co2_apriori_levels|20>'
end KEYWORDS

begin FILE
   dest_filename = <co2_apriori_filename>

   begin LOAD_CO2_APRIORI
      sounding_info_file = <case_soundinginfo_filename>
      surface_type       = <brdf_type>
      global_mean_xco2   = '<global_mean_xco2|>'
      num_levels         = '<num_co2_apriori_levels>'
   end LOAD_CO2_APRIORI
end FILE

begin FILE
   source_filename = oco_l2.inp

   begin PICK 
      section  = PARAMETER_DEFINITION->GAS 
      keyword  = a_priori
      which_section = CO2 
      template = <co2_apriori_filename>
   end PICK 

end FILE