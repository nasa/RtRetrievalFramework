PRO barplot2, x, y, colors=colors, xoffset=xoffset, _extra=_extra

; Like BARPLOT, but lets you have many different "categories" of y.
; stacks the categories on top of one another
;
; INPUT VARIABLES
; x : the x locations of the bars (vector)
; y : the y locations of the bars (*,n_types)

; KEYWORDS
; WIDTH : the width of the bars, as a fraction of the normalized x width.
; XOFFSET: the x-offset for each bar, in the same units as WIDTH.
; OPLOT : Overplot on the existing plot (you may want to set the OFFSET keyword for this).
; _EXTRA: Pretty much any keyword allowed to PLOT and OPLOT (such as COLOR, XTITLE, YTITLE, etc).

	ny = n_elements(y[0,*])
	if n_elements(colors) ne ny then colors = findgen(ny)/(ny-1.)*254
	if n_elements(xoffset) eq 0 then xoffset=0.
	barplot, x, y[*,0], _extra=_extra, xoff=xoffset
	barplot, x, y[*,0], _extra=_extra, /oplot, col=colors[0], xoff=xoffset

	ytot = y[*,0] ; cumulative amount of y
	for i =1, n_elements(y[0,*])-1 do begin
		barplot, x,ytot + y[*,i],yoffset=ytot, _extra=_extra, $
			xoff=xoffset, /oplot, col = colors[i]
		ytot += y[*,i]
	endfor


END