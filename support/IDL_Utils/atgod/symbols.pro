;+
; <pre>
; NAME:	    
;           SYMBOLS
;
; PURPOSE:
;	    Create custom plotting symbols
;
; CATEGORY:
;           Plotting
;
; CALLING SEQUENCE:
;	    symbols,symbol_number,scale	
; INPUTS:
;	    symbol_number=
;                1 = open circle
;                2 = filled circle
;                3 = arrow pointing right
;                4 = arrow pointing left
;                5 = arrow pointing up
;                6 = arrow pointing down
;                7 = arrow pointing up and left (45 degrees)
;		 8 = arrow pointing down and left
;		 9 = arrow pointing down and right.
;		 10 = arrow pointing up and right.
;		 11 through 18 are bold versions of 3 through 10
;		 19 = horizontal line
;		 20 = box
;		 21 = diamond
;		 22 = triangle
;		 30 = filled box
;		 31 = filled diamond
;		 32 = filled triangle
;		 
;	    Scale = size of symbols.
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
;		Color = color of symbols
;               Thick = thickness of lines drawn to construct symbol.
; OUTPUTS:
; OPTIONAL OUTPUT PARAMETERS:
; COMMON BLOCKS:
; SIDE EFFECTS:
;           The desired symbol is stored in the user buffer and 
;	    will be plotted if !P.PSYM = 8.
;
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;		Jeff Bennett, U of Colorado, 198?
;            Added to idlmeteo from the Windt-Library   Nov-1992  oet\@sma.ch
;
; modified by MLS - 3/12/96
;
; </pre>
;-
pro symbols,nsym,scale,color=col,thick=thick
on_error,2
fill = 0
thick = 3
case 1 of
     (nsym le 2):   begin                         ;circles
                      ;for large scales increase number of points for res.
                      if scale ge 4 then a = findgen(25) else $
                      a = findgen(13)
                      a = a * (3.14159 / 6.)       ;(0 - 12 or 24) pi/6
                      xarr = cos(a)
                      yarr = sin(a)
                      if nsym eq 2 then fill = 1
                    end
     ((nsym ge 3)*(nsym le 18)):   begin           ;arrow heads
                      xarr = fltarr(5)
                      yarr = xarr
                      xarr(1) = 10.
                      xarr(2) = 6.
                      yarr(2) = 2.
                      ;nsyms greater than 10 should be filled arrows
                      if nsym gt 10 then begin
                         xarr(3) = 6. 
                         xarr(4) = 10.
                         yarr(3) = -2.
                         fill = 1
                      endif else begin
                         xarr(3) = 10.
                         xarr(4) = 6.
                         yarr(4) = -2.
                      endelse
                      case 1 of
                         (nsym eq 3): dummy = 0b
                         (nsym eq 4): xarr = -1.*xarr
                         ((nsym eq 11)+(nsym eq 12)): begin
                            xarr = extrac(xarr,0,11)
                            yarr = extrac(yarr,0,11)
                            yarr(6) = 0.5
                            xarr(7) = 6
                            yarr(7) = 0.5
                            xarr(8) = 6
                            yarr(8) = -0.5
                            yarr(9) = -0.5
                            if nsym eq 12 then begin
                               rotation,xarr,yarr,180,nx,ny
                               xarr = nx
                               yarr = ny
                            endif
                                                     end
                         ((nsym eq 5)+(nsym eq 13)): begin
                            temp = xarr
                            xarr = yarr
                            yarr = temp
                                                     end
                         ((nsym eq 6)+(nsym eq 14)): begin
                            temp = -1.*xarr
                            xarr = yarr
                            yarr = temp
                                                     end
                         ((nsym ge 7)*(nsym le 10) + $
                         (nsym ge 15)*(nsym le 18)): begin
                            case 1 of
                                 ((nsym eq 7)+(nsym eq 15)): deg = 45
                                 ((nsym eq 8)+(nsym eq 16)): deg = 135
                                 ((nsym eq 9)+(nsym eq 17)): deg = 225
                                 ((nsym eq 10)+(nsym eq 18)): deg = 315
                            endcase
                            rotation,xarr,yarr,deg,nx,ny
                            xarr = nx
                            yarr = ny
                                                  end   ;end nsym ge 7
                      endcase
                                   end    ;nsym between 3 and 18
     ((nsym eq 20)+(nsym eq 21)+(nsym eq 30)+(nsym eq 31)):  begin
                      xarr = fltarr(5) + 3
                      yarr = xarr
                      xarr(1) = -3.
                      xarr(2) = -3.
                      yarr(2) = -3.
                      yarr(3) = -3.
                      if (nsym eq 21)+(nsym eq 31) then begin
                         rotation,xarr,yarr,45,nx,ny
                         nx = 0.70 * nx     ;shrink the x direction
                         xarr = nx
                         yarr = ny
                      endif
                      if nsym ge 30 then fill = 1
                                   end    ;nsym 20,21,30,31
     ((nsym eq 22)+(nsym eq 32)):  begin  ;side length 6, 0 at centroid
                      yarr = fltarr(4) - 6./4.
                      xarr = fltarr(4) - 6./2.
                      xarr(1) = 6./2.
                      xarr(2) = 0.
                      yarr(2) = 6.*sqrt(3.)/2. - 6./4.
                      if nsym eq 32 then fill = 1
                                    end
     else:                          begin
                      xarr = fltarr(2) + 1
                      yarr = xarr * 0.
                      xarr(1) = -1.
                                    end
endcase
;
xarr = xarr * scale
yarr = yarr * scale
;
;set symbol buffer
if keyword_set(col) then usersym,xarr,yarr,fill=fill,color=col,thick=thick $
   else usersym,xarr,yarr,fill=fill,thick=thick
;
return
end
