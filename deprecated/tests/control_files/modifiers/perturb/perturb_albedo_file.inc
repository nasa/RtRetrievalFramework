begin FILE
  source_filename = <case_albedo_filename>
  only_if  = 'eval( "<brdf_type>".lower() == "lambertian" )'

  ########
  # ALBEDO

  begin MODIFY_DATA
      columns = 0
      modify  = 'abs( <original> * 0.85 )'
  end MODIFY_DATA

  begin MODIFY_DATA
      columns = 1
      modify  = 'abs( <original> * 1.05 )'
  end MODIFY_DATA

  begin MODIFY_DATA
      columns = 2
      modify  = 'abs( <original> * 0.80 )'
  end MODIFY_DATA

end FILE
