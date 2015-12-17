; $Id: allocatecmykcolor.pro,v 1.4 2008/01/31 00:06:04 livesey Exp $

function AllocateCMYKColor, c, m, y, k, percent=percent
common CMYKMappings, mappings

factor = 1.0
if keyword_set(percent) then factor = 0.01

cc = factor * ( c + k )
mm = factor * ( m + k )
yy = factor * ( y + k )
mx = max ( [ cc, mm, yy ] )
if mx gt 1.0 then begin
  cc = cc / mx
  mm = mm / mx
  yy = yy / mx
endif

r = 1.0 - cc
g = 1.0 - mm
b = 1.0 - yy

result = AllocateColor ( 255*r, 255*g, 255*b )

;; Now record this color for later
thisMapping = { TCMYKMapping, $
  rgb: [r,g,b], cmyk:[c,m,y,k]*factor }
if n_tags(mappings) eq 0 then mappings = thisMapping $
else mappings = [ mappings, thisMapping ]

return, result

end
