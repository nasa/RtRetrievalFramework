pro FlushCMYKMappings
common CMYKMappings, mappings
oneMapping = { TCMYKMapping, $
  rgb: fltarr(3), cmyk:fltarr(4) }
mappings = replicate ( oneMapping, 2 )
;; Do black
mappings(0).rgb  = 0.0
mappings(0).cmyk = [ 0.0, 0.0, 0.0, 1.0 ]
;; Do white
mappings(1).rgb  = 1.0
mappings(1).cmyk = 0.0

end
