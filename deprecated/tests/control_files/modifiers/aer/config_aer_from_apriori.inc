include = <CONTROL_FILE_DIRECTORY>/standard_aerosol_types.inc

begin FILE
  source_filename = oco_l2.inp

  #########
  # AEROSOL

  begin CREATE_AEROSOL_BLOCKS
      use_buffered_file = True
      from_file  = <case_aerosol_filename>
      set_retrieval_vector = true
      profile_names = '<use_aerosol_profile_names|>'
      type_names = '<use_aerosol_type_names|>'
  end CREATE_AEROSOL_BLOCKS

  # Need to fix mie/mom file names since the standard profile names
  # do not have matching mie/mom files for them
  # The following only matches when standard aerosol types are setup
  # from the aerosol apriori file
  begin SUBSTITUTE
     from = '\/wc\.m'
     to   = '/<wc_type>.m'
     keyword = mie_file moment_file
  end SUBSTITUTE

  begin SUBSTITUTE
     from = '\/ic\.m'
     to   = '/<ic_type>.m'
     keyword = mie_file moment_file
  end SUBSTITUTE

  begin SUBSTITUTE
     from = '\/aero1\.m'
     to   = '/<aero1_type>.m'
     keyword = mie_file moment_file
  end SUBSTITUTE

  begin SUBSTITUTE
     from = '\/aero2\.m'
     to   = '/<aero2_type>.m'
     keyword = mie_file moment_file
  end SUBSTITUTE

end FILE


