;+
; NAME:
;        CSYM
;
; PURPOSE:
;        This function defines a standard sequence of plotting symbols
;        that can be used in place of PSYM e.g. in a call to
;        PLOT. Example: plot,x,y,psym=csym(1) draws filled circles.
;        A chart of these new symbols can be created with Showsym.
;
; CATEGORY:
;        General Graphics
;
; CALLING SEQUENCE:
;        PLOT,X,Y,PSYM=SYM(NUMBER)
;
; INPUTS:
;        NUMBER    ->   symbol number
;
;		see code below
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;        function returns the symbol number to be used with PSYM= in the
;        PLOT command
;
; SUBROUTINES:
;        SHOWCSYM : Can be used to produce a symbol chart for reference
;        (Type .r sym, then showsym, optionally with the /PS option).
;        Extra keywords are passed to PLOTS, so you can e.g. choose
;        a fancy color for your chart.
;
; REQUIREMENTS:
;
; NOTES:
;        This function produces a side effect in that the USERSYM procedure
;        is used to create a symbol definition. It's meant for usage within
;        the PLOT, OPLOT, etc. command
;
;
; MODIFICATION HISTORY:
;        mgs, 22 Aug 1997: VERSION 1.00
;        mgs, 10 Sep 1999: - added SHOWSYM procedure
;        mgs, 26 Aug 2000: - changed copyright to open source
;
;-
;

;###########################################################################


pro showcsym,ps=ps,_EXTRA=e, thick=thick

   FORWARD_FUNCTION CSYM


   psflag = keyword_set(PS)
   if (psflag) then begin
      olddev = !D.NAME
      set_plot,'PS'
      device,/COLOR,bits=8,xsize=8,ysize=5,yoffset=3,/INCHES, $
           filename='symbols.ps'
   endif

   c = csym(allowed = a)
   na = n_elements(a)
   plot,findgen(na),/NODATA,xstyle=4,YSTYLE=4
   for i=0,na-1 do begin
      plots,1,na-i,PSYM=CSYM(a[i], thick=thick),_EXTRA=e
      xyouts,0.5,na-i-0.2,strtrim(a[i],2),align=1.
   endfor

   if (psflag) then begin
      device,/close
      set_plot,olddev
      print,'Symbollist created as symbols.ps.'
   endif

   return
end

; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function csym,number, allowed=allowed, thick=thick

     on_error,2  ; return to caller


	 allowed = [indgen(15), 16,17,18,19]
     if(n_elements(number) eq 0) then return,0  ; default

     result=8    ; default: return psym=8, i.e. user defined symbol


; define some help variables for
; circle :
     phi=findgen(48)*(!PI*2/48.)
     phi = [ phi, phi(0) ]


     case abs(number) of

         0  : result = 0    ; no symbol
         1	: usersym, [ -1, 1, 0, 0, 0], [ 0, 0, 0, 1, -1 ], thick=thick	; plus
         2	: usersym,  1.2*[0, 0, 0,0,-.7,.7,0,.7,-.7], 1.2*[0, 1, -1,0,-.7,.7,0,-.7,.7 ], thick=thick	; star
         3	: result = 3	; dot
		 4  : usersym, [ 0, 1, 0, -1, 0 ], [ 1, 0, -1, 0, 1 ], thick=thick
                            ; open diamond

		 5	: usersym, [ -1, 0, 1, -1 ], [ -1, 1, -1, -1 ], thick=thick
                            ; open upward triangle

		 6	: usersym, [ -1, 1, 1, -1, -1 ], [ 1, 1, -1, -1, 1 ], thick=thick
                            ; open square

		 7	: result = 7	; X

		 8  : usersym, [ -1, 0, 1, -1 ], [  1, -1, 1, 1 ], thick=thick
                           ; open downward triangle
         9 	: usersym, cos(phi), sin(phi), thick=thick
                            ; open circle
         10 : usersym, [ -1, 1, -1, -1 ], [1, 0, -1, 1 ], thick=thick
                           ; rightfacing triangle, open

         11  : usersym, [ 1, -1, 1, 1 ], [1, 0, -1, 1 ], thick=thick
                           ; leftfacing triangle, open

		 12	: usersym, [ 0, 1, 0, -1, 0 ], [ 1, 0, -1, 0, 1 ], /fill, thick=thick
                            ; filled diamond

 		 13 : usersym, [ -1, 0, 1, -1 ], [ -1, 1, -1, -1 ], /fill, thick=thick
                            ; filled upward triangle

		 14 : usersym, [ -1, 1, 1, -1, -1 ], [ 1, 1, -1, -1, 1 ], /fill, thick=thick
                            ; filled square

		 15 : usersym, [ -0.5, 0.5, 0.5, -0.5, -0.5], [-1,-1,1,1,-1], /fill, thick=thick
         16 : usersym, [ -1, 0, 1, -1 ], [  1, -1, 1, 1 ], /fill, thick=thick
                           ; filled downward triangle
         17	: usersym, cos(phi), sin(phi), /fill, thick=thick
                            ; filled circle

         18 : usersym, [ -1, 1, -1, -1 ], [1, 0, -1, 1 ], /fill, thick=thick
                           ; rightfacing triangle, filled

         19 : usersym, [ 1, -1, 1, 1 ], [1, 0, -1, 1 ], /fill, thick=thick
                           ; leftfacing triangle, filled


       else : begin
              print,'invalid symbol number - set to 0'
              result = 0
              end

     endcase

	if number LT 0 then result = -1 * result
return,result
end

