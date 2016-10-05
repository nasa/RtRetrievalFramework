; $Id: setps.pro,v 1.14 2005/07/14 01:06:43 livesey Exp $
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; This routine setups all the various things associated with ps
; output.
;
; @keyword portrait     {type=Boolean} {Default=0}
;                       Sets up page in portrait mode.
; @keyword landscape    {type=Boolean} {Default=0}
;                       Sets up page in landscape mode.  Is ignored if
;                       portrait is also set.  Causes unpredicatble
;                       results when encapsulated is set.
; @keyword full         {type=Boolean} {Default=0}
;                       Full portrait page.
; @keyword big          {type=Boolean} {Default=0}
;                       11x17 page.  Is ignored when full is set.
; @keyword encapsulated {type=Boolean} {Default=0}
;                       Sets up page in encapsulated postscript.
; @keyword color        {type=Boolean} {Default=0}
;                       Allows page to be written in color.
; @keyword download     Does not do anything.  Left in for legacy reasons.
; @keyword margin       Does not do anything.  Left in for legacy reasons.
; @keyword caslon       {type=Boolean} {Default=0}
;                       Sets font as acaslon-Regular.
; @keyword optima       {type=Boolean} {Default=0}
;                       Sets font as Optima.
; @keyword gillsans     {type=Boolean} {Default=0}
;                       Sets font as GillSans
; @keyword filename     {type=String} {Default=idl.ps}
;                       The name of the file being selected
; @keyword noThick      {type=Boolean} {Default=0}
;                       When unset, all of the plot thick settings are
;                       at 2 instead of 1. 
; @keyword hOffset      {type=Float} {Default=1.5}
;                       The size of the margins on the left and right
;                       hand of the postscript page in centimeters.
;                       Only used when either full or big are set.
;                       This takes into account the orientation of
;                       the page. 
; @keyword vOffset      {type=Float} {Default=1.5}
;                       The size of the margins on the top and bottom
;                       of the postscript page in centimeters.
;                       Only used when either full or big are set.
;                       This takes into account the orientation of
;                       the page.                
;
; @author Nathaniel Livesey; January 18, 1999
; @version $Revision: 1.14 $ $Date: 2005/07/14 01:06:43 $
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro SetPs, portrait=portrait, landscape=landscape, full=full, big=big, $
           encapsulated=encapsulated, color=color, download=download, $
           margin=margin, caslon=caslon, optima=optima, gillSans=gillSans, $
           filename=filename, noThick=noThick, hOffset=hOffset, vOffset=vOffset

Set_Plot, 'ps'
Device, /CLOSE, encapsulated=Keyword_Set(encapsulated), $
        color=color, bits=4+4*Keyword_Set(color)
IF Keyword_Set(filename) THEN Device, filename=filename

CASE 1 OF
  Keyword_Set(optima) : Device,user_font='Optima'
  Keyword_Set(caslon) : Device,user_font='ACaslon-Regular'
  Keyword_Set(gillsans) : Device,user_font='GillSans'
  ELSE : Device, /TIMES
ENDCASE
!P.font=0
IF NOT Keyword_Set(noThick) THEN BEGIN
  !P.thick = 2.0
  !X.thick = 2.0
  !Y.thick = 2.0
ENDIF ELSE BEGIN
  !P.thick = 1.0
  !X.thick = 1.0
  !Y.thick = 1.0
ENDELSE

; Now do the page size type stuff

;margin=keyword_set(margin) ; unused
landscape = Keyword_Set(landscape) AND NOT Keyword_Set(portrait)
portrait = landscape EQ 0
big = Keyword_Set(big)
full= Keyword_Set(full)


IF full EQ 1 OR big EQ 1 THEN BEGIN

  hOffset = N_Elements(hOffset) EQ 1 ? hOffset[0] : 1.5
  vOffset = N_Elements(vOffset) EQ 1 ? vOffset[0] : 1.5

  ;; The width and height will be in cm
  IF full EQ 1 THEN BEGIN
    pageWidth = 8.5 * 2.54
    pageHeight = 11.0 * 2.54
  ENDIF ELSE BEGIN              ; big EQ 1
    pageWidth = 11.0 * 2.54
    pageHeight = 17.0 * 2.54
  ENDELSE

  IF portrait EQ 1 THEN BEGIN
    xSize = pageWidth - 2.0 * hOffset
    ySize = pageHeight - 2.0 * vOffset
    xOffset = hOffset
    yOffset = vOffset
  ENDIF ELSE BEGIN              ; landscape
    xSize = pageHeight - 2.0 * hOffset
    ySize = pageWidth - 2.0 * vOffset
    xOffset = vOffset
    yOffset = hOffset + xSize
  ENDELSE
ENDIF

;; [x|y]Size and [x|y]OffsetAdj are all undefined when neither big or
;; full are set.
Device, xSize=xSize, ySize=ySize, xOffset=xOffset, yOffset=yOffset, $
        landscape=landscape, portrait=portrait

END

; $Log: setps.pro,v $
; Revision 1.14  2005/07/14 01:06:43  livesey
; Added AdobeCalson
;
; Revision 1.13  2004/07/06 23:30:06  fullerr
; Rewrote to make simpler, added in hOffset and vOffset keywords
;
; Revision 1.12  2004/01/02 18:24:26  livesey
; Added GillSans options
;
; Revision 1.11  2003/09/18 23:47:46  fullerr
; Reformatted/cleaned up code, and added documentation
;
; Revision 1.10  2002/08/27 19:31:37  livesey
; Better handling of noThick
;
; Revision 1.9  2002/07/30 22:14:09  fullerr
; Hopefully fixed bug with switching between portrait and landscape
;
; Revision 1.8  2002/06/10 22:27:56  livesey
; Made portrait default
;
; Revision 1.7  2002/03/12 23:07:21  livesey
; Various changes
;
; Revision 1.6  2001/07/03 00:08:48  livesey
; Reformatted and indented etc.
;
; Revision 1.5  2001/01/23 05:25:32  livesey
; Sent home

; Revision 1.4  2000/10/03 16:07:50  livesey
; Removed atmqty

; Revision 1.3  2000/06/22 01:40:58  livesey
; Regular commit

; Revision 1.2  1999/10/20 00:41:14  livesey
; Added font stuff

; Revision 1.1  1999/01/18 20:21:28  livesey
; Initial revision


