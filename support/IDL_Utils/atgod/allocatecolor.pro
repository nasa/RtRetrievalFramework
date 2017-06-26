; $Id: allocatecolor.pro,v 1.7 2003/11/13 23:31:08 fullerr Exp $
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; This routine allocates a single color requested by the user, it
; returns the value for this color.
; <p>
; Part of the ColorBoss suite.
;
; @param red          {type=Long} {Required}
;                     The red component of the desired color.
; @param green        {type=Long} {Required}
;                     The green component of the desired color.
; @param blue         {type=Long} {Required}
;                     The blue component of the desired color.
; @keyword noTransfer {type=Boolean} {Default=0}
;                     Explicitly disallows using the transfer function
;                     on the allocated color.
; @keyword rgb        {type=Long}
;                     The color in a single value.  Overrides the
;                     values of red, green, and blue.
; @keyword hex        {type=Boolean} {Default=0}
;                     Does the color in hex.
;
; @returns The value of the newly allocated color.
;
; @author Nathaniel Livesey
; @version $Revision: 1.7 $ $Date: 2003/11/13 23:31:08 $
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function AllocateColor, red, green, blue, $
                        noTransfer=noTransfer, rgb=rgb, hex=hex
common ColBoss,config

IF N_Elements(config) EQ 0 THEN InitColorBoss

IF N_Elements(rgb) ne 0 THEN BEGIN
  IF keyword_set(hex) THEN BEGIN
    red = rgb MOD 256L
    green = rgb / 256L 
    blue = green / 256L
    green = green MOD 256L
  ENDIF ELSE BEGIN
    blue = rgb MOD 256L
    green = rgb / 256L
    red = green / 256L
    green = green MOD 256L
  ENDELSE
ENDIF

; First check that this is not straight black or white, if it is we'll
; ignore the transfer function
ignoreTransfer = 0

IF (red EQ 0) and (green EQ 0) and (blue EQ 0) THEN ignoreTransfer = 1
IF (red EQ 255) and (green EQ 255) and (blue EQ 255) THEN ignoreTransfer = 1

useRed = red
useGreen = green
useBlue = blue
IF config.transferFunction NE '' AND NOT Keyword_Set(noTransfer) AND $ 
  ignoreTransfer EQ 0 THEN BEGIN
  dummy = Call_Function(config.transferFunction,$
                        red=useRed,green=useGreen,blue=useBlue)
ENDIF

IF config.pseudo THEN BEGIN
  IF config.firstFree EQ !D.n_colors THEN Message, 'No colors left'

  ;; Sort out result, update information
  result = config.firstFree
  config.firstFree = config.firstFree + 1

  ;; Read table in, set this value, write table back
  Tvlct, tableRed, tableGreen, tableBlue, /GET

  tableRed[result] = useRed
  tableGreen[result] = useGreen
  tableBlue[result] = useBlue

  TVLCT, tableRed, tableGreen, tableBlue

ENDIF ELSE BEGIN

  ;; In true or direct color, then this is easy
  result = Long(useRed) + '100'xl * Long(useGreen) + $
           '10000'xl * Long(useBlue)
ENDELSE

Return, result
END

; $Log: allocatecolor.pro,v $
; Revision 1.7  2003/11/13 23:31:08  fullerr
; Added documentation and did some reformatting
;
; Revision 1.6  2001/10/18 00:50:57  livesey
; Added rgb and hex options
;
; Revision 1.5  2001/07/03 00:08:48  livesey
; Reformatted and indented etc.
;
; Revision 1.4  1999/10/22 01:09:18  livesey
; Nightly commit

; Revision 1.3  1999/10/21 17:18:39  livesey
; Added transfer function stuff to colorboss

; Revision 1.2  1999/01/18 18:35:43  livesey
; Americanised the spellings.

; Revision 1.1  1999/01/18 18:32:48  livesey
; Initial revision



