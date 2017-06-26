begin FILE
  source_filename = <case_aerosol_filename>

  begin MERGE_PROFILES
     wc    = 'wc_.*'
     ic    = 'ic_.*'
     aero1 = 'cont'    'clat_0' 'clat_3' 'clat_4' 'clat_5'
     aero2 = 'oceanic' 'clat_1' 'clat_2' 'clat_6'
  end MERGE_PROFILES


end FILE
