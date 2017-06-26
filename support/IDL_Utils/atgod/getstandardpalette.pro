; $Id: getstandardpalette.pro,v 1.14 2008/07/29 18:00:34 livesey Exp $
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; This is a useful routine to get a standard palette set.  Set the
; desired colors as boolean keywords, ie GetStandardPalette, /BLACK.
; <p>
; The valid colors are:
; <ul>
; <li>Black
; <li>White
; <li>Gray (grey)
; <li>LightGrey (lightgray)
; <li>LightLightGrey (lightlightgray)
; <li>DarkGrey (darkgray)
; <li>Red
; <li>Green
; <li>Blue
; <li>Cyan
; <li>Magenta
; <li>Yellow
; <li>Orange
; </ul>
; <p>
; Part of the ColorBoss suite
;
; @keyword set    {type=Boolean} {Default=0}
;                 When given, sets the color and background color.
; @keyword paper  {type=Boolean} {Default=0}
;                 When set, use color condusive to paper (ie black
;                 ink, white background.  The default is reversed).
;                 Only works when the keyword set is set.
; @keyword _extra A list of keywords that are the colors to allocate.
;                 Their values don't matter.
;
; @returns The desired palette set.
;
; @author Nathaniel Livesey
; @version $Revision: 1.14 $ $Date: 2008/07/29 18:00:34 $
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function GetStandardPalette, set=set, paper=paper, _extra=requests, cmyk=cmyk, all=all

set = Keyword_Set(set)
nNames = N_Tags(requests)
IF nNames EQ 0 THEN BEGIN
  names = 'DUMMY'
  nNames = 1
ENDIF ELSE names = Tag_Names(requests)

IF set EQ 1 THEN BEGIN  ; Make sure black and white are on the list.
  IF (where(names EQ 'BLACK'))(0) EQ -1 THEN BEGIN
    names = [names,'BLACK']
    nNames = nNames + 1
  ENDIF
  IF (where(names EQ 'WHITE'))(0) EQ -1 THEN BEGIN
    names = [names,'WHITE']
    nNames = nNames + 1
  ENDIF
ENDIF

FOR tag = 0, nNames - 1 DO BEGIN
  found = 1
  IF KEYWORD_SET ( cmyk ) THEN BEGIN
    CASE names[tag] OF
      'BLACK'          : color = [0,  0,  0,  100]
      'WHITE'          : color = [0,  0,  0,  0  ]

      'GREY'           : color = [0,  0,  0,  50 ]
      'LIGHTGREY'      : color = [0,  0,  0,  25 ]
      'LIGHTLIGHTGREY' : color = [0,  0,  0,  12 ]
      'DARKGREY'       : color = [0,  0,  0,  75 ]
      'ALMOSTWHITE'    : color = [0,  0,  0,  1  ]
      
      'GRAY'           : color = [0,  0,  0,  50 ]
      'LIGHTGRAY'      : color = [0,  0,  0,  25 ]
      'LIGHTLIGHTGRAY' : color = [0,  0,  0,  12 ]
      'DARKGRAY'       : color = [0,  0,  0,  75 ]
      
      'RED'            : color = [0,  100,80 ,0  ]
      'GREEN'          : color = [100,0  ,100,10 ]
      'BLUE'           : color = [100,60 ,0,  10 ]

      'CYAN'           : color = [100,0,  0,  0  ]
      'MAGENTA'        : color = [0,  100,0,  0  ]
      'YELLOW'         : color = [0,  0  ,100,0  ]

      'LIGHTRED'       : color = [0,  50 ,50 ,0  ]
      'LIGHTGREEN'     : color = [50 ,0,  50 ,0  ]
      'LIGHTBLUE'      : color = [50 ,25 ,0,  0  ]

      'DARKRED'        : color = [0,  100,70 ,40 ]
      'DARKGREEN'      : color = [100,0,  80 ,40 ]
      'DARKBLUE'       : color = [100,80, 0,  70 ]

      'LIGHTCYAN'      : color = [50 ,0,  0,  0  ]
      'LIGHTMAGENTA'   : color = [0,  50 ,0,  0  ]
      'LIGHTYELLOW'    : color = [0,  0,  50 ,0  ]

      'DARKCYAN'       : color = [100,0,  0,  50  ]
      'DARKMAGENTA'    : color = [0,  100,0,  50  ]
      'DARKYELLOW'     : color = [0,  0,  100,50  ]

      'LIGHTLIGHTRED'  : color = [0,  25, 25, 0  ]
      'LIGHTLIGHTGREEN': color = [25, 0,  25, 0  ]
      'LIGHTLIGHTBLUE' : color = [25, 12, 0,  0  ]

      'ORANGE'         : color = [0,  50 ,100,0  ]
      'BROWN'          : color = [15, 75 ,100,30 ]
      'PURPLE'         : color = [100,100,0,  0  ]
      'DUMMY'          : found = 0
      ELSE             : found = 0
    ENDCASE
    IF found EQ 1 THEN $
      thisColor = AllocateCMYKColor ( color[0], color[1], color[2], color[3], /percent )
  ENDIF ELSE BEGIN
    CASE names[tag] OF
      'BLACK'          : color = [0  ,0  ,0  ]
      'WHITE'          : color = [255,255,255]
      
      'GREY'           : color = Replicate(128, 3)
      'LIGHTGREY'      : color = Replicate(192, 3)
      'LIGHTLIGHTGREY' : color = Replicate(224, 3)
      'DARKGREY'       : color = Replicate(64, 3)
      'ALMOSTWHITE'    : color = Replicate(254, 3)
      
      'GRAY '          : color = Replicate(128, 3)
      'LIGHTGRAY'      : color = Replicate(192, 3)
      'LIGHTLIGHTGRAY' : color = Replicate(224, 3)
      'DARKGRAY'       : color = Replicate(64, 3)
      
      'RED'            : color = [255,0  ,0  ]
      'GREEN'          : color = [0  ,255,0  ]
      'BLUE'           : color = [0  ,0  ,255]

      'CYAN'           : color = [0  ,255,255]
      'MAGENTA'        : color = [255,0  ,255]
      'YELLOW'         : color = [255,255,0  ]

      'LIGHTRED'       : color = [255,128,128]
      'LIGHTGREEN'     : color = [128,255,128]
      'LIGHTBLUE'      : color = [128,128,255]

      'DARKRED'        : color = [128,0  ,0  ]
      'DARKGREEN'      : color = [0  ,128,0  ]
      'DARKBLUE'       : color = [0  ,0,  128]

      'LIGHTCYAN'      : color = [128,255,255]
      'LIGHTMAGENTA'   : color = [255,128,255]
      'LIGHTYELLOW'    : color = [255,255,128]

      'DARKCYAN'       : color = [0  ,128,128]
      'DARKMAGENTA'    : color = [128,0  ,128]
      'DARKYELLOW'     : color = [128,128,0  ]

      'LIGHTLIGHTRED'  : color = [255,192,192]
      'LIGHTLIGHTGREEN': color = [192,255,192]
      'LIGHTLIGHTBLUE' : color = [192,192,255]

      'ORANGE'         : color = [255,128,0  ]
      'BROWN'          : color = [120,60, 0  ]
      'PURPLE'         : color = [74, 36, 94 ]

      'DUMMY'          : found = 0
      ELSE             : found = 0
    ENDCASE
    IF found EQ 1 THEN $
      thisColor = AllocateColor(color[0], color[1], color[2])
  ENDELSE
  IF found EQ 1 THEN BEGIN
    result = N_Elements(result) EQ 0 ? $
      Create_Struct(names[tag], thisColor) : $
      Create_Struct(result, names[tag], thisColor)
  ENDIF ELSE IF names[tag] NE 'DUMMY' THEN BEGIN
    Message, 'Unable to find color '+ StrLowCase(names[tag])
  ENDIF
ENDFOR

IF set EQ 1 THEN BEGIN
  IF Keyword_Set(paper) THEN BEGIN
    !P.color = result.black
    !P.background = result.white
  ENDIF ELSE BEGIN
    !P.color = result.white
    !P.background = result.black
  ENDELSE
END

Return, result
END

; $Log: getstandardpalette.pro,v $
; Revision 1.14  2008/07/29 18:00:34  livesey
; Added AlmostWhite (for fills that don't get deleted by
; epstopdf/pdflatex)
;
; Revision 1.13  2008/03/02 20:05:59  livesey
; Added more colors etc.
;
; Revision 1.12  2007/01/12 00:44:50  livesey
; Added more colors
;
; Revision 1.11  2006/08/27 18:54:30  livesey
; Added lightmagenta
;
; Revision 1.10  2005/01/07 01:48:56  livesey
; Added some more colors
;
; Revision 1.9  2004/11/09 15:22:34  livesey
; Added light colors
;
; Revision 1.8  2003/11/14 17:24:14  fullerr
; Bug Fix
;
; Revision 1.7  2003/11/14 01:34:57  fullerr
; Added documentation and did some reformatting
;
; Revision 1.6  2003/11/14 01:34:26  fullerr
; Added documentation and did some reformatting
;
; Revision 1.5  2002/06/03 00:43:36  livesey
; Sent home
;
; Revision 1.4  2001/09/02 05:17:37  livesey
; Regular commit
;
; Revision 1.3  2001/08/07 14:36:10  livesey
; Sent stuff home
;
; Revision 1.2  2001/07/03 00:08:48  livesey
; Reformatted and indented etc.
;
; Revision 1.1  1999/01/18 18:35:56  livesey
; Initial revision


