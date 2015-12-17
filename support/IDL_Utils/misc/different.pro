function different, x, z=z, index=index, old=old, sorted=sorted

; given an array x of values, this function returns a (shortened) array of the
; x-values that are different from each other.
;

; 2 Different Methods:
;
; 1. NEW METHOD
;		Uses "sort" trick of craig Marquardt.  Pros : Can work on any sortable array (not just integers).
;		Just as fast (or faster) then OLD METHOD.  Cons : Cannot tell you how many repeats each unique element had
;	(although could probably be modified to do so).

; 2.  OLD METHOD
;     The values of x MUST be integers. operates on the baglean-sort principle.
;
; KEYWORDS
;		OLD : Set this to use the "baglean-sort" method.
; 	z : an array of the # of different entries of each existing element of the returned array (/OLD only)
;		INDEX : Set this to return the indices of the unique elements (NEW), or for "z" to equal these indices (OLD).

if keyword_set(old) then begin
	x0 = min(x)
	x1 = max(x)

	N = (x1-x0) + 1

	z = lonarr(N)
	if ~keyword_set(index) then begin
		for i = 0L, n_elements(x)-1 do begin
			argi = x[i] - x0
			z[argi] = z[argi] + 1
		endfor
		wz = where(z) ; elements of z with some hits
		z = z[wz] ; # of different entries of each existing x
		return, wz + x0
	endif else begin
		argi = long(x[0] - x0)
		z[argi] = 1
		w = 0
		for i = 1L, n_elements(x)-1 do begin
			argi = x[i] - x0
			z[argi] = z[argi] + 1
			if z[argi] eq 1 then w = [w,i]
		endfor
		 z = w
	   return, x[z]
	endelse
endif else begin
	; Craig Marquardt's Method:
	if n_elements(x) LE 1 then if keyword_set(index) then return, 0L else return, x[0]
 
  if keyword_set(sorted) then begin
  	ii = lindgen(n_elements(x))
	b = temporary(x)
  endif else begin
    ii = sort(x) & b = x(ii) ; b is the sorted version of x
  endelse
  wh = where(b NE shift(b, +1L), ct) ; find unique values of b
  if keyword_set(sorted) then x = temporary(b)
  if ct GT 0 then if keyword_set(index) then return, ii[wh] else return, x[ii[wh]]

  if keyword_set(index) then return, 0L else return, x[0]

endelse
end
