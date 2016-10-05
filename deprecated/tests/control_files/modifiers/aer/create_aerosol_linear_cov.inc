begin FILE
   source_filename = <case_aerosol_filename>
   dest_filename = <aerosol_covariance_filename>

   begin CREATE_SIMPLE_COV_MODULE
      columns = '1'
      scaling = '0.01'
   end CREATE_SIMPLE_COV_MODULE
end FILE

begin FILE
  source_filename = oco_l2.inp

  begin PICK
      section = PARAMETER_DEFINITION->AEROSOL
      keyword = covariance
      template = <aerosol_covariance_filename>
  end PICK

end FILE
