pro Plot_Generic_Residuals, x, y1, y2, NAMES=names, TITLE=plot_title, XTITLE=title, YTITLE=ytitle

   ;; Set up for plot parts
   if keyword_seT(names) then $
     leg_names = names[0:1] $
   else $
     leg_names = ['y1', 'y2']

   if not keyword_set(plot_title) then $
     plot_title = ''

   if not keyword_set(x_title) then $
     x_title = ''

   if not keyword_set(y_title) then $
     y_title = ''

   txt_colors = [1, 1]
   item_colors = [1, 2]
   item_psym = [0, 0]
   item_lines = [0, 2]

   charsize_orig = !p.charsize
   ;!p.charsize = 0.5

   legend_spacing = 0
               
   multi_orig = !p.multi
   x_style_orig = !x.style
   y_style_orig = !y.style

   !p.multi = [0, 0, 2, 0, 0]
   !x.style = 1
   !y.style = 1  


   residual = 100. * ((y2 - y1) / max(y1))
   smax = max([max(y1), max(y2)])
   rmax = max([abs(min(residual)), max(residual)])

   plot, x, y1, /nodata, color=1, position=[0.13, 0.10, 0.95, 0.82], $
         xrange = [min(x), max(x)], $
         xtickname =['', '', '', '', '', '', '', '', '', ''], yrange = [0., smax], $
         xtitle=xtitle, ytitle=ytitle
   oplot, x, y1, linestyle=item_lines[0], psym=item_psym[0], color = item_colors[0]
   oplot, x, y2, linestyle=item_lines[1], psym=item_psym[1], color = item_colors[1]
   
   Legend_DR, leg_names, TEXTCOLORS=txt_colors, COLORS=item_colors, PSYM=item_psym, $
              LINESTYLE=item_lines, CORNERS=corners, POS=[2,2], SPACING=legend_spacing
   
   xydims = [corners[2]-corners[0],corners[3]-corners[1]]
   chdim=[!d.x_ch_size/float(!d.x_size),!d.y_ch_size/float(!d.y_size)]
   pos = [!x.window[1]-chdim[0]/2.0-xydims[0], !y.window[0]+xydims[1]/2.0+chdim[1]]
   
   polyx  = [ pos[0], pos[0] + xydims[0], pos[0] + xydims[0], pos[0] ]
   polyy  = [ pos[1], pos[1], pos[1] - xydims[1], pos[1] - xydims[1] ]
   
   ;; Fill the area where the legend will be drawn so that no data will
   ;; leak through and obscure the legend
   polyfill, polyx, polyy, /NORM, COLOR=0
   
   Legend_DR, leg_names, TEXTCOLORS=txt_colors, COLORS=item_colors, PSYM=item_psym, $
              LINESTYLE=item_lines, POS=pos, SPACING=legend_spacing
   

   plot, x, residual, /nodata, color = 1, title=plot_title, $
         position = [0.13, 0.83, 0.95, 0.95], xtitle = ' ', $
         xtickname =[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], $
         ytitle = 'Residual (%)', yrange = [-1*rmax, rmax]
   
   oplot, x, dblarr(n_elements(residual)), linestyle = 0, color = 1 ; zero line
   oplot, x, residual, linestyle = 0, color = 1
   
   !p.multi = multi_orig
   !x.style = x_style_orig
   !y.style = y_style_orig

   !p.charsize = charsize_orig

end

