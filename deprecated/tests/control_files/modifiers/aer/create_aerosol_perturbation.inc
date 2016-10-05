# Aerosol must be in linear mode for finite difference jacobians
begin FILE
   source_filename = <case_aerosol_filename>
   dest_filename   = <aerosol_perturb_filename>

   begin MODIFY_DATA
      columns     = '1, 999'
      modify      = '1e-10'
   end MODIFY_DATA

end FILE

begin FILE
  source_filename = oco_l2.inp

  begin PICK
      section = PARAMETER_DEFINITION->AEROSOL
      keyword = perturb
      template = <aerosol_perturb_filename>
  end PICK

end FILE
