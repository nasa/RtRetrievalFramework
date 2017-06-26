; $Id: array_nextdivision.pro,v 1.7 2004/11/09 15:22:05 livesey Exp $
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; This routine moves onto the next division on the screen (if there is
; one!)
; <p>
; Part of the Array_plots suite.
; This set of idl routines makes a better job of the !P.multi type
; things that idl can do.  It allows the user to divide and
; recursively subdivide the screen/page into sections.
;
; @author Nathaniel Livesey
; @version $Revision: 1.7 $ $Date: 2004/11/09 15:22:05 $
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro array_nextDivision
common arrayPlots,setup

N = N_Elements(setup)-1
setup[N].current = setup[N].current + 1
IF setup[N].current EQ setup[N].no THEN BEGIN ; The last plot in section
  IF N GT 0 THEN BEGIN          ; This has real parents
    setup=setup[0:N-1]
    Array_NextDivision          ; Move onto the next pane in the parent.
  ENDIF ELSE BEGIN              ; This is the ultimate parent, i.e. whole page
    setup[0].current = 0
    !P.region = setup[0].currentRegion
    !P.position = 0
  ENDELSE
ENDIF ELSE BEGIN                ; Not the last plot, move onto the next one.
  index = setup[N].current
  setup[N].currentRegion = [setup[N].x0[index], $
                            setup[N].y0[index], $
                            setup[N].x0[index] + setup[N].xSize[index], $
                            setup[N].y0[index] + setup[N].ySize[index]]
  IF setup[N].showEdges THEN BEGIN
    Plots, /Normal, setup[N].currentRegion([0,2,2,0,0]), $
           setup[N].currentRegion([1,1,3,3,1])
  END
  IF setup[N].byPosition THEN BEGIN
    !P.position = setup[N].currentRegion
    !P.region = 0
  ENDIF ELSE BEGIN
    !P.region = setup[N].currentRegion
    !P.position = 0
  ENDELSE
ENDELSE

END

; $Log: array_nextdivision.pro,v $
; Revision 1.7  2004/11/09 15:22:05  livesey
; Bug fix, better handling of end game
;
; Revision 1.6  2003/09/20 00:22:23  fullerr
; Reformatted/cleaned up code, and added documentation
;
; Revision 1.5  2001/07/03 00:08:48  livesey
; Reformatted and indented etc.
;
; Revision 1.4  2000/06/16 02:14:02  livesey
; Regular commit

; Revision 1.3  2000/05/19 22:19:22  livesey
; Reset whichever of !P.position and !P.region we're not using each time.

; Revision 1.2  2000/05/19 20:35:29  livesey
; Dealt with the new byPosition functionality

; Revision 1.1  1999/01/18 19:13:28  livesey
; Initial revision


