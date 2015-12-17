begin FILE
  source_filename = <case_aerosol_filename>

  # Make sure first the data is not already in log (most likely < 0.0 if so
  begin MODIFY_DATA
     columns = '1,99'
     modify = '<original> >= 0.0 and math.log( <original> < 1e-20 and 1e-20 or <original> ) or <original>'
     not_if = '<use_original_aerosol_mode|False>'
   end MODIFY_DATA

  begin PICK
     section  = HEADER
     keyword  = retrieval_mode
     template = logarithmic
     not_if = '<use_original_aerosol_mode|False>'
  end PICK

end FILE
