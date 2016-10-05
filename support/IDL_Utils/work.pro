PRO Work, $
  COLORTABLE=colorTable

  ;; Expand these dirs into the IDL search path
  dirsToExpand = [ getenv('L2_SUPPORT_PATH') + '/IDL_Utils', $
                   '/groups/gds/deliv/Idl/Container/latest/pro', $
                   '/groups/gds/deliv/Idl/Bitfields/latest/pro', $
                   '/groups/gds/deliv/Idl/FillValues/latest/pro', $
                   '/groups/gds/deliv/Idl/OcoHdf/latest/pro', $
                   '/groups/gds/deliv/Idl/IdlUnit/latest/unittest/pro' $
                 ]

  ;; Select the path delimiter character.
  CASE !D.NAME OF
      'VMS' : delim = ','
      'WIN' : delim = ';'
      ELSE    : delim = ':'
  ENDCASE

  ;; Expand and added to !PATH all directories with .pro files.
  FOR iPath=0, N_Elements(dirsToExpand)-1L DO $
            !PATH= Expand_Path('+' + dirsToExpand[iPath]) + delim + !PATH

  Set_Plot, 'X'

  ;; Initialize the device for color display.
  IF ( (!D.N_Colors GE 256L) AND $
       ((!D.Name EQ 'X') OR (!D.NAME EQ 'MacOS')) ) THEN $
    Device, Pseudo_Color=8

  Device, Decomposed=0, Bypass_Translation=1, Retain=2

  ;; Save the current color table.
  TVLCT, savedR, savedG, savedB, /GET
  colorTable = [[savedR],[savedG],[savedB]]
  
  ;; Load the base color table.
  LoadCT, 22

  ;; Load the working colors
  Setup_Color_Table

END ;; End of work
