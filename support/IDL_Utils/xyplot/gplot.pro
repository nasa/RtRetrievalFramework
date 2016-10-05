pro gplot,xx,ex,yy,ey,index,xcp,ycp,symfile,captions,nsym,ipanel,npanel,orientation,xtxt,ytxt,text,ntxt
;
; Procedure to plot arrays XX, YY with their respective error bars EX, EY
; and with different symbol/colors/linestyles according to array INDEX.
;
; Inputs:
;      XX(*)        Array of x-values
;      EX(*)        Array of x-uncertainties
;      YY(*)        Array of y-values
;      EY(*)        Array of y-uncertainties
;      INDEX(MSYM)  Array of indices at which the symbol/color is changed
;      XCP(*)       Array of x-ordinates of captions
;      YCP(*)       Array of y-ordinates of captions
;      SYMFILE      Name of file (e.g xyplot.sym) containing list of symbols
;      CAPTIONS     String array containing caption text
;      NSYM         Number of different symbols to be used.
;      IPANEL       Index of the current plot panel.
;      NPANEL       Total number of plot panels.
;      ORIENTATION  Orientation of plot (either "portrait" or "landscape")
;      XTXT(NTXT)   Array of x-ordinates of text strings
;      YTXT(NTXT)   Array of y-ordinates of text strings
;      TEXT(NTXT)   String array of text strings
;      NTXT         Number of text strings to be plotted
;
; Outputs:
;      none
;
;
;  Notes:  Text strings are printed in black
;          Captions are plotted alongside their symbol in the same color
;
ccc=string(80)
charsize=1.0
captsize=1.0
linethick=1.0
colortable=0
nrow=1
ncol=1
orientation=string("landscape")
openr,unit,symfile,/get_lun  ;  Open file of symbols.
readf,unit,nhead,nfields
for k=2,nhead do begin
   readf,unit,ccc
   if strpos(ccc,"multi_panel:") eq 0 then reads,strmid(ccc,12,99),ncol,nrow
   if strpos(ccc,"orientation:") eq 0 then reads,strmid(ccc,13,99),orientation
   if strpos(ccc,"color_table:") eq 0 then reads,strmid(ccc,12,99),colortable
   if strpos(ccc,"char_size:") eq 0   then reads,strmid(ccc,10,99),charsize
   if strpos(ccc,"capt_size:") eq 0   then reads,strmid(ccc,10,99),captsize
   if strpos(ccc,"line_thick:") eq 0  then reads,strmid(ccc,11,99),linethick
endfor
npanel=nrow*ncol
!p.thick=linethick
!p.charsize=charsize/sqrt(npanel)
captsize=captsize/sqrt(npanel)
if colortable ge 0 then loadct,colortable ; -ve value forces black-white plots
!p.multi=0
if npanel gt 1  then !p.multi=[(npanel-ipanel) mod npanel,ncol,nrow]
;
; Draw frame and plot first point
plot, xx(0:0), yy(0:0)
;map_set,/cylindrical,/continent,latdel=5.,londel=5,/grid,limit=[60,10,75,45]
;map_continents,/countries,/hires
;map_set,/cylindrical,/continent,latdel=5.,londel=5,/grid,limit=[50,-135,70,-75]
;map_continents,/countries,/hires
;
; Write filename and time in lower left corner
xyouts,0,0,text(0),alignment=0.0,charsize=0.4*captsize,/normal
;
; Write text strings to the appropriate locations on the plot.
for itxt=1,ntxt-1 do  xyouts,xtxt(itxt),ytxt(itxt),text(itxt),charsize=captsize
;
; Correct y-offset between symbols and caption text
dy=nrow*(!y.range(1)-!y.range(0))*captsize/200
klo=index(0)
for jsym=0,nsym-1 do begin
  khi=index(jsym+1)-1
; Read symbol parameters and define its array
  if not eof(unit) then repeat readf,unit,ccc until (eof(unit) or strmid(ccc,0,1) ne ";")
  if not eof(unit) then reads, ccc, size,ntimes,nls,nrot,fill,col,ls,xeb,yeb
  size=captsize*size
  a = 2*!pi*(ntimes*findgen(nls+1)/nls+0.125*nrot)
  if fill eq 1 then usersym,size*sin(a),size*cos(a),/fill else usersym,size*sin(a),size*cos(a)
;
  if(khi ge klo) then begin
    !p.psym = 0
    if (colortable ge 0) then !p.color=col
    for i=klo,khi  do begin                  ;  plot error bars (individually)
    if(xeb eq 1) then oplot,[xx(i)-ex(i),xx(i)+ex(i)],[yy(i),yy(i)] ;x-error bar
    if(yeb eq 1) then oplot,[xx(i),xx(i)],[yy(i)-ey(i),yy(i)+ey(i)] ;y-error bar
    if(xeb ge 2) then oplot,[xx(i)-ex(i),xx(i)+ex(i)],[yy(i)-ey(i),yy(i)+ey(i)]
    endfor
    if(ls ge 0) then oplot,xx(klo:khi),yy(klo:khi),linestyle=ls  ;connecting line
    if(size gt 0) then oplot,xx(klo:khi),yy(klo:khi), psym=8     ;  plot data points
    if(strlen(captions(jsym)) gt 0) then begin  ; draw symbols & captions
       oplot,xcp(jsym:jsym),ycp(jsym:jsym),psym=8
       xyouts, xcp(jsym), ycp(jsym)-dy,'   '+captions(jsym),charsize=captsize
    endif
  endif
  klo=khi+1
endfor ; jsym=0,nsym-1
;
close, unit
free_lun, unit

return
end
