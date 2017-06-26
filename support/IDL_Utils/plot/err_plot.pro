PRO ERR_PLOT, X, YLOW, YHIGH, WIDTH=WIDTH, FLIP=flip, $
  _EXTRA=EXTRA_KEYWORDS

;- Check arguments
if (n_params() ne 3) then message, $
  'Usage: ERR_PLOT, X, YLOW, YHIGH'
if (n_elements(x) eq 0) then $
  message, 'Argument X is undefined'
if (n_elements(ylow) eq 0) then $
  message, 'Argument YLOW is undefined'
if (n_elements(yhigh) eq 0) then $
  message, 'Argument YHIGH is undefined'

;- Check keywords
if (n_elements(width) eq 0) then width = 0.02

;- Plot the error bars
for index = 0L, n_elements(x) - 1L do begin

  ;- Plot vertical bar using data coordinates
  if keyword_set(flip) then begin
      ydata = [x[index], x[index]]
      xdata = [ylow[index], yhigh[index]]
  endif else begin
      xdata = [x[index], x[index]]
      ydata = [ylow[index], yhigh[index]]
  endelse

  plots, xdata, ydata, /data, noclip=0, $
         _extra=extra_keywords

  ;- Compute horizontal bar width in normal coordinates
  normalwidth = (!x.window[1] - !x.window[0]) * width

  ;- Plot horizontal bar using normal coordinates

  if keyword_set(flip) then begin
      lower = convert_coord(ylow[index], x[index], $
                            /data, /to_normal)
      upper = convert_coord(yhigh[index], x[index], $
                            /data, /to_normal)

      xdata1 = [lower[0], lower[0]]
      xdata2 = [upper[0], upper[0]]
      ylower = [lower[1] - 0.5 * width, lower[1] + 0.5 * width]
      yupper = [upper[1] - 0.5 * width, upper[1] + 0.5 * width]
  endif else begin
      lower = convert_coord(x[index], ylow[index], $
                            /data, /to_normal)
      upper = convert_coord(x[index], yhigh[index], $
                            /data, /to_normal)

      xdata1 = [lower[0] - 0.5 * width, lower[0] + 0.5 * width]
      xdata2 = [lower[0] - 0.5 * width, lower[0] + 0.5 * width]
      ylower = [lower[1], lower[1]]
      yupper = [upper[1], upper[1]]
  endelse

  plots, xdata1, ylower, /normal, noclip=0, $
         _extra=extra_keywords
  plots, xdata2, yupper, /normal, noclip=0, $
         _extra=extra_keywords


endfor

END
