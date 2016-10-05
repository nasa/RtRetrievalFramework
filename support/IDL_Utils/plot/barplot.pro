PRO barplot, x, y, width=width, xoffset=xoffset, yoffset=yoffset, oplot=oplot, _extra=_extra

; Makes a bar plot.  This is similar to the IDL routine BAR_PLOT, except that it is more like PLOT;
;  it has an x axis, and the bars are placed appropriately on the x axis.   Accepts keywords like xlog
;  and ylog, and pretty much any other keywords accepted by PLOT and OPLOT.

; INPUT VARIABLES
; x : the x locations of the bars
; y : the y locations of the bars

; KEYWORDS
; WIDTH : the width of the bars, as a fraction of the normalized x width.
; OFFSET: the x-offset for each bar, in the same units as WIDTH.
; OPLOT : Overplot on the existing plot (you may want to set the OFFSET keyword for this).
; _EXTRA: Pretty much any keyword allowed to PLOT and OPLOT (such as COLOR, XTITLE, YTITLE, etc).

	N = n_elements(x)
	if n_elements(width) eq 0 then width = 1./(2.*N)
	oplot = keyword_set(oplot)
	if n_elements(xoffset) eq 0 then begin
		if oplot then xoffset = 2*width else xoffset = 0.0
	endif

	case n_elements(yoffset) of
	0 : yoffsets = fltarr(n)
	1 : yoffsets = fltarr(n) + yoffset
	n : yoffsets = yoffset
	endcase

	; make the initial plot if necessary
	if ~oplot then plot, x, y, /nodata, _extra=_extra

	xlog = !x.type
	ylog = !y.type

	if xlog then x_ = alog10(x) else x_ = x

	y0 = !y.crange[0] ; lowest value on y axis.

	w = width * (!x.crange[1] - !x.crange[0])
	o = xoffset* (!x.crange[1] - !x.crange[0])

	for i = 0, N-1 do begin
    	y0i = y0 + yoffsets[i]
		if ylog then y0i = 10.^y0i
		yi = y[i] > y0i
		xi = x_[i]
		yvec = [y0i,y0i,yi,yi]
		xvec = [xi-w, xi+w, xi+w, xi-w] + o
		if xlog then xvec = 10.^xvec
		if y[i] GT y0i then polyfill, xvec, yvec, _extra=_extra
	endfor

END