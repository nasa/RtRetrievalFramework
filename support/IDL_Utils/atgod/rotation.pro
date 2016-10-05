;+
;
; return (in x/yout) the coordinates of the input points after
; rotation through 'angle' degrees
; 
; Meant to replace the missing routine `rotation', which is
; called by 'symbols', which is from some library which seems to have
; disappeared. I'm only guessing what this routine did, but I'm tired
; of having symbols fail because of the lack of this routine.
; 
; @parameter xin {input}{type=float vector}{required}
;   The x coordinates of the points before rotation
; @parameter yin {input}{type=float vector}{required}
;   The y coordinates of the points before rotation
; @parameter angle {input}{type=float vector}{required}
;   The angle to rotate the points. 
; @parameter xout {output}{type=float vector}{required}
;   The x coordinates of the points after rotation
; @parameter yout {output}{type=float vector}{required}
;   The y coordinates of the points after rotation
;
; @keyword radians {type=boolean}{optional}{default=false}
;   assume 'angle' is in degrees unless /radians
; 
;-
PRO rotation, xin,yin,angle,xout,yout,radians=radians
   
   IF n_params() LT 5 THEN BEGIN 
     Message,'Need all 5 arguments!',/info
     return
   ENDIF 
   tangle = keyword_set(radians) EQ 0? angle/!radeg : angle
   xout = cos(tangle)*xin-sin(tangle)*yin
   yout = sin(tangle)*xin + cos(tangle)*yin
return
END

; $Id: rotation.pro,v 1.1 2007/09/24 21:30:23 whdaffer Exp $
;
; Mod Log :
; $Log: rotation.pro,v $
; Revision 1.1  2007/09/24 21:30:23  whdaffer
; Initial (and, I hope only) revision
;
;


