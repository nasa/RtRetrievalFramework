; $Id: array_newdivision.pro,v 1.16 2006/04/18 16:30:59 fullerr Exp $
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; This routine is called to set up a new division of the screen,
; it can be a completely new screen or a sub division of a parent
; division.
; <p>
; Part of the Array_plots suite.
; <p>
; I should probably change this some time to use pointers, rather than
; fixed array sizes. Then again, if it's not broken....
;
; @author Nathaniel Livesey
; @version $Revision: 1.16 $ $Date: 2006/04/18 16:30:59 $
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; @hidden
; A subprogram of Array_NewDivision, takes care of everything but
; outputting the title. 
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO Array_NewDivision_Driver, nx=nx, ny=ny, byPosition=byPosition, $
                              xRatio=xRatio, yRatio=yRatio, $
                              xMarginRatio=xMarginRatio, $ 
                              yMarginRatio=yMarginRatio, $
                              sense=sense, noErase=noErase, $
                              showEdges=showEdges, reset=reset, ignore=ignore

  COMPILE_OPT IDL2, HIDDEN
  COMMON arrayPlots, setup

  MAX_DIVISIONS = 50L
  MAX_PLOTS = MAX_DIVISIONS ^ 2

  IF N_Elements(setup) EQ 0 OR Keyword_Set(reset) THEN BEGIN
    ;; Initialise the system
    setup={ nx:1,ny:1,no:1,$
            current:0,$
            $ ;; Used to mark areas we should not look at for charsize
            ignore:0B, $
            parentRegion:[0.0,0.0,1.0,1.0],$
            $ ;; This addition comes in handy to figure out depth
            parentRegionString:'', $ 
            currentRegion:[0.0,0.0,1.0,1.0], $
            byPosition:0,$
            showEdges:0, $
            xSize:FltArr(MAX_PLOTS), $
            ySize:FltArr(MAX_PLOTS), $
            x0:FltArr(MAX_PLOTS), $
            y0:FltArr(MAX_PLOTS) $
          }
    setup.parentRegionString = StrJoin(StrTrim(setup.parentRegion, 2), ',')
  ENDIF

  depth = N_Elements(setup)
  
  IF N_Elements(nx) EQ 0 THEN nx=1
  IF N_Elements(ny) EQ 0 THEN ny=1
  IF N_Elements(xRatio) EQ 0 THEN xRatio = FltArr(nx)+(1.0/nx)
  IF N_Elements(yRatio) EQ 0 THEN yRatio = FltArr(ny)+(1.0/ny)
  IF N_Elements(xMarginRatio) EQ 0 THEN xMarginRatio = FltArr(nx+1)
  IF N_Elements(yMarginRatio) EQ 0 THEN yMarginRatio = FltArr(ny+1)
  ;; By default sense is 1 unless nx==1
  IF N_Elements(sense) EQ 0 THEN sense = nx EQ 1 ? 2 : 1

  nXMarginRatio = N_Elements(xMarginRatio)
  nYMarginRatio = N_Elements(yMarginRatio)
  IF nXMarginRatio EQ 1 THEN BEGIN
    nXMarginRatio = nx + 1
    xMarginRatio = Replicate(xMarginRatio[0], nXMarginRatio)
  ENDIF
  IF nYMarginRatio EQ 1 THEN BEGIN
    nYMarginRatio = ny + 1
    yMarginRatio = Replicate(yMarginRatio[0], nYMarginRatio)
  ENDIF
  ;; Now do sanity checks on the calling code.
  CASE 1 OF
    nx GT MAX_DIVISIONS      : Message, 'Too many X divisions'
    ny GT MAX_DIVISIONS      : Message, 'Too many Y divisions'
    Min(xRatio) LT 0.0       : Message, 'Inappropriate xRatio'
    Min(yRatio) LT 0.0       : Message, 'Inappropriate yRatio'
    Min(xMarginRatio) LT 0.0 : Message, 'Inappropriate xMarginRatio'
    Min(yMarginRatio) LT 0.0 : Message, 'Inappropriate yMarginRatio'
    nXMarginRatio NE nx+1    : Message, 'Inappropriate size for xMarginRatio'
    nYMarginRatio NE ny+1    : Message, 'Inappropriate size for yMarginRatio'
    ELSE:                       ; Do nothing
  ENDCASE

  IF depth EQ 1 AND NOT Keyword_Set(noErase) THEN Erase

  ;; Now deal with the title lines (moved up to caller)
  ;; noLines=N_Elements(title)

  ;; IF N_Elements(myCharsize) EQ 1 AND noLines NE 0 THEN $
  ;;   myCharsize = Replicate(myCharsize, noLines)
  ;; maxChar = Max(myCharsize)

  ;;titleHeight= $
  ;;   Convert_Coord([0, maxchar* !D.y_ch_size * noLines], /DEVICE, /TO_NORMAL)
  ;; titleHeight=titleHeight[1]

  ;; region = setup[depth-1].currentRegion

  ;; IF N_Elements(alignment) EQ 1 AND noLines NE 0 THEN $
  ;;   alignment = Replicate(alignment, noLines)

  ;; yTitle=region[3] - titleHeight
  ;; FOR line=0,noLines-1 DO BEGIN
  ;;   CASE alignment[noLines-line-1] OF
  ;;     0.5: xTitle=(region[2] + region[0]) / 2.0
  ;;     0.0: xTitle = region[0]
  ;;     1.0: xTitle = region[2]
  ;;   ENDCASE
  ;; ENDFOR

  ;;   thisLine = title[noLines-line-1]
  ;;   yLine = yTitle+line*titleHeight/noLines
  ;;   thisSize = myCharsize[noLines-line-1]
  ;;   IF StrPos(thisLine,'##') EQ -1 THEN BEGIN
  ;;     XYOuts, xTitle, yLine, /NORMAL, $
  ;;       thisLine, alignment=alignment[noLines-line-1], charsize=thisSize
  ;;   ENDIF ELSE BEGIN
  ;;     titleSections = StrSplit(thisLine, '##', /EXTRACT)
  ;;     XYOuts, region[0], yLine, /NORMAL, $
  ;;             titleSections(0), alignment=0.0, charsize=thisSize
  ;;     XYOuts, (region[0]+region[2])/2.0, yLine, /NORMAL, $
  ;;             titleSections[1], alignment=0.5, charsize=thisSize
  ;;     XYOuts, region[2], yLine, /normal, $
  ;;             titleSections[2], alignment=1.0, charsize=thisSize
  ;;   ENDELSE
  ;; ENDFOR
  ;; titleHeight = 0

  setup = [setup, setup[depth-1]]
  setup[depth].nx = nx
  setup[depth].ny = ny
  setup[depth].current = -1
  setup[depth].parentRegion = setup[depth-1].currentRegion ; - [0.0,0.0,0.0,titleHeight]
  setup[depth].parentRegionString = StrJoin(StrTrim(setup[depth].parentRegion, 2), ',')
  setup[depth].currentRegion = [0.0,0.0,0.0,0.0]
  setup[depth].ignore = Keyword_Set(ignore)
  setup[depth].byPosition = Keyword_Set(byPosition)
  setup[depth].showEdges = Keyword_Set(showEdges) OR (setup[(depth-1)>0].showEdges)
  setup[depth].no = nx * ny

  xFullSize = setup[depth].parentRegion[2] - setup[depth].parentRegion[0]
  yFullSize = setup[depth].parentRegion[3] - setup[depth].parentRegion[1] ; - titleHeight

  totXRatio = Total(xRatio) + Total(xMarginRatio)
  totYRatio = Total(yRatio) + Total(yMarginRatio)
  xSize = xFullSize * xRatio / TotXRatio
  ySize = yFullSize * yRatio / TotYRatio
  xMarginSize = xFullSize * xMarginRatio / TotXRatio
  yMarginSize = yFullSize * yMarginRatio / TotYRatio
  
  xCumul = FltArr(nx*2+1)
  FOR i = 0, nx-1 DO BEGIN
    xCumul[2*i]   = xMarginSize[i]
    xCumul[2*i+1] = xSize[i]
  ENDFOR
  xCumul[2*nx] = xMarginSize[nx]

  yCumul=FltArr(ny*2+2)
  yCumul[0] = setup[depth].parentRegion[3] ; - titleHeight
  FOR i=0,ny-1 DO BEGIN
    yCumul[2*i+1] = -yMarginSize[i]
    yCumul[2*i+2] = -ySize[i]
  ENDFOR

  xCumul = Total(xCumul, /CUMULATIVE)
  yCumul = Total(yCumul, /CUMULATIVE)

  xStart = xCumul[IndGen(nx)*2]
  yStart = yCumul[IndGen(ny)*2+2]


  ;; Now sort out the ordering of the plots, select row or col. major
  inds = IndGen(setup[depth].no)
  IF sense EQ 1 THEN BEGIN
    setup[depth].x0 = xStart[inds MOD nx] + setup[depth].parentRegion[0]
    setup[depth].y0 = yStart[inds / nx] ;+setup[depth].parentRegion[1]
    setup[depth].xSize = xSize[inds MOD nx]
    setup[depth].ySize = ySize[inds / nx]
  ENDIF ELSE BEGIN
    setup[depth].x0 = xStart[inds / ny] + setup[depth].parentRegion[0]
    setup[depth].y0 = yStart[inds MOD ny] ;+setup[depth].parentRegion[1]
    setup[depth].xSize = xSize[inds / ny]
    setup[depth].ySize = ySize[inds MOD ny]
  ENDELSE
  IF Max(Abs(setup[depth].y0)) GT 1.0 THEN stop

  ;; Now call the NextDivision routine to set up the screen.

  Array_NextDivision

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; This routine is called to set up a new division of the screen,
; it can be a completely new screen or a sub division of a parent
; division.
; <p>
; Part of the Array_plots suite.
; <p>
; I should probably change this some time to use pointers, rather than
; fixed array sizes. Then again, if it's not broken....
;
; @author Nathaniel Livesey
; @version $Revision: 1.16 $ $Date: 2006/04/18 16:30:59 $
;
; @keyword title        {type=String} {Default=''}
;                       Title string (or string array for multi-lines)
; @keyword charsize     {type=Float} {!P.charsize}
;                        Size for title
; @keyword alignment    {type=Float}
;                       Alignment for title (as in xyouts)
; @keyword byPosition   {type=Boolean} {default=0}  
;                       Divide using !p.position, not !p.region
; @keyword nx           {type=Int}
;                       Number of divisions in the x direction.
; @keyword ny           {type=Int}
;                       Number of divisions in the x direction.
; @keyword xRatio       {type=Float[]} 
;                       Relative sizes of lengthwise divisions (equal assumed)
; @keyword yRatio       {type=Float[]} 
;                       Relative sizes of height divisions (equal assumed)
; @keyword xMarginRatio {type=Float[]} {Default=0}
;                       Relative sizes of xMargins (nx+1)
; @keyword yMarginRatio {type=Float[]} {Default=0}
;                       Relative sizes of yMargins (ny+1)
; @keyword sense        {type=Int} {Default=1}
;                       Go clockwise or anticlockwise (1 or 2)
; @keyword noErase      {type=Boolean} {default=0}  
;                       Don't erase, even on reset
; @keyword showEdges    {type=Boolean} {default=0}  
;                       Diagnostic, draw box to show edges
; @keyword reset        {type=Boolean} {default=0}  
;                       Start with fresh screen, even if child expctd
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro Array_NewDivision, title=title, charsize=charsize, alignment=alignment, $
                       byPosition=byPosition, nx=nx, ny=ny, $
                       xRatio=xRatio, yRatio=yRatio, $
                       xMarginRatio=xMarginRatio, yMarginRatio=yMarginRatio, $
                       sense=sense, noErase=noErase, $
                       showEdges=showEdges, reset=reset
  COMPILE_OPT IDL2
  COMMON arrayPlots, setup

  IF Keyword_Set(reset) EQ 1 THEN BEGIN
    Array_NewDivision_Driver, /RESET, noErase=noErase, /IGNORE
  ENDIF

  ;; First figure out the title
  depth = N_Elements(setup)
  noLines = N_Elements(title)

  IF noLines GT 0 THEN BEGIN    ; There is a title
    ;; Figure out charsize/alignment
    nCharSize = N_Elements(charSize)
    IF nCharSize EQ 0 THEN BEGIN
      IF depth GT 0 THEN BEGIN
        inds = Where(setup.ignore EQ 0, cnt)
        temp = (cnt EQ 0 ? 0 : $
                N_Elements(Uniq(setup[inds].parentRegionString, $
                                Sort(setup[inds].parentregionstring)))-1)
      ENDIF ELSE temp = 0
      myCharsize = Replicate((1.5-0.5*temp)>0.5, noLines)
    ENDIF ELSE IF nCharSize NE noLines THEN BEGIN
      myCharSize = Replicate(charSize[0], noLines)
    ENDIF ELSE BEGIN
      myCharSize = charSize
    ENDELSE

    hasAlign = 1
    nAlign = N_Elements(alignment)
    IF nAlign EQ 0 THEN BEGIN
      hasAlign = 0
      myAlign = Replicate(0.5, noLines)
    ENDIF ELSE IF nAlign NE noLines THEN BEGIN
      myAlign = Replicate(alignment[0], noLines)
    ENDIF ELSE BEGIN
      myAlign = alignment
    ENDELSE

    ;; Now figure out how this fits into the big picture of the plot.  A
    ;; single 1.0 charsize character is estimated to take up about 1/60th
    ;; of the page.  This is conservative and says each normal line takes
    ;; up 1/50th of the page, use this to figure out the ratio of title
    ;; versus plot.
    titlePortion = Total(myCharSize)
    plotPortion = 50 - titlePortion

    ;; Now divide up the page plot  
    Array_NewDivision_Driver, nx=1, ny=2, yRatio=[titlePortion, plotPortion], $
                              /BYPOSITION, /IGNORE

    ;; Now put in the plot into the first part
    Array_NewDivision_Driver, ny=noLines, yRatio=myCharSize, /BYPOSITION, /IGNORE

    ;; Now start plotting in the title line by line
    FOR i = 0, noLines - 1 DO BEGIN
      ;; Check to see for columns
      IF StrPos(title[i], '##') EQ -1 THEN BEGIN ; 1 column
        Plot, [0], xRange=[0,1], yRange=[0,1], xStyle=5, yStyle=5, $
              /NOERASE, /NODATA
        XYOuts, myAlign[i], 0, title[i], charSize=myCharSize[i], align=myAlign[i]
        Array_NextDivision
      ENDIF ELSE BEGIN
        ;; Divide columns
        cols = StrSplit(title[i], '##', /REGEX, /EXTRACT, count=nCols)
        ;; If nCols == nx then preserve margin ratios
        matchMargins = 0
        IF N_Elements(xMarginRatio) GT 0 THEN BEGIN
          If nCols EQ nx THEN BEGIN
            xmr = xMarginRatio
            matchMargins = 1
          ENDIF ELSE xmr = 0
        ENDIF ELSE xmr = 0

        Array_NewDivision_Driver, nx=nCols, xMarginRatio=xmr, /BYPOSITION, /IGNORE
        FOR j = 0, nCols - 1 DO BEGIN
          ;; Figure out alignment
          IF NOT matchMargins THEN BEGIN
            CASE j OF
              0: thisAlign = 0.0 
              nCols - 1: thisAlign = 1.0
              ELSE: thisAlign = myAlign[i]
            ENDCASE
          ENDIF ELSE thisAlign = myAlign[i]
          
          Plot, [0], xRange=[0,1], yRange=[0,1], xStyle=5, yStyle=5, $
                /NOERASE, /NODATA
          XYOuts, thisAlign, 0, cols[j], charSize=myCharSize[i], align=thisAlign
          Array_NextDivision
        ENDFOR
      ENDELSE
    ENDFOR
  ENDIF

  ;; Now do the rest of the business
  Array_NewDivision_Driver, byPosition=byPosition, nx=nx, ny=ny, $
                            xRatio=xRatio, yRatio=yRatio, $
                            xMarginRatio=xMarginRatio, yMarginRatio=yMarginRatio, $
                            sense=sense, noErase=noErase, showEdges=showEdges

END
  
; $Log: array_newdivision.pro,v $
; Revision 1.16  2006/04/18 16:30:59  fullerr
; Bug Fix
;
; Revision 1.15  2006/04/12 23:03:34  fullerr
; Fixed bug in automatic charsize
;
; Revision 1.14  2003/09/29 16:13:50  fullerr
; Bug Fix
;
; Revision 1.13  2003/09/22 20:52:21  fullerr
; Cleaned up code, did some reformatting and documentation
;
; Revision 1.12  2003/03/13 22:07:13  fullerr
; A new and (hopefully) much improved method for inserting titles
;
; Revision 1.11  2002/03/12 23:07:21  livesey
; Various changes
;
; Revision 1.10  2002/01/10 01:55:19  livesey
; Regular commit
;
; Revision 1.9  2001/11/03 01:18:02  livesey
; Regular commit
;
; Revision 1.8  2001/07/04 02:14:15  livesey
; Increased limits
;
; Revision 1.7  2001/07/03 00:08:48  livesey
; Reformatted and indented etc.
;
; Revision 1.6  2001/07/02 23:15:05  livesey
; Regular commit

; Revision 1.5  2000/06/16 02:13:58  livesey
; Regular commit

; Revision 1.4  2000/06/14 02:29:18  livesey
; Various changes notably to atgod_hdf

; Revision 1.3  2000/05/19 22:18:57  livesey
; Renormalised default xRatio and yRatio

; Revision 1.2  2000/05/19 20:35:14  livesey
; Added the margin and byPosition functionality

; Revision 1.1  1999/01/18 19:12:19  livesey
; Initial revision



