pro standard_testplots, summary_file, stats_file, out_filename, log_file, $
  plot_title=plotTitle, $
  max_chi2=max_chi2, max_aod_tr=max_aod_tr, $
  max_aod_rt=max_aod_rt, max_xco2_sig=max_xco2_sig, $
  max_rat12=max_rat12, max_rat32=max_rat32

   if not Keyword_Set(max_chi2) then max_chi2=0.5
   if not Keyword_Set(max_aod_tr) then max_aod_tr=0.3
   if not Keyword_Set(max_aod_rt) then max_aod_rt=0.3
   if not Keyword_Set(max_xco2_sig) then max_xco2_sig=1.
   if not Keyword_Set(max_rat12) then max_rat12=0.5
   if not Keyword_Set(max_rat32) then max_rat32=0.5

   if not Keyword_Set(plotTitle) then $
     plotTitle='L2 Full Physics Run Summary'

   plotTitle = plotTitle + ': '
   
   ; read in relevant info from various output files
   ; concatenate good and bad points together
   print, 'Loading summary file for plotting: ', summary_file
   dataG   = read_fp_summary_dat(summary_file)
   print, 'Loading stats file for plotting: ', stats_file
   statsG  = read_fp_summary_dat(stats_file)

   ncases=n_elements(dataG)

   ;; Catch errors with setting up data values. Simply return with
   ;; warning, because certain required contents of the input files was
   ;; likely missing
   Catch, error_status
   if error_status ne 0 then begin
       print, 'An error occured when setting up data values for plotting:'
       print, !ERROR_STATE.MSG
       print, out_filename, ' will not be created.'
       Catch, /CANCEL
       return
   endif
   
   if Keyword_Set(log_file) then begin
       print, 'Loading orbit simulator log file for plotting: ', log_file
       logdata = read_orbit_simulator_log(log_file)

       where_fp_id = where(tag_names(logdata) eq "FOOTPRINT_ID")

       ;; Make sure we can find the row we need before looking for it
       if where_fp_id[0] ne -1 then begin
           ;; find the cases in log file which were retrieved and had good_quality
           ;; put their index values in idx_log
           idx_log=intarr(ncases)
           ids=ulong64(dataG.case_number)
           for i=0,ncases-1 do begin
               ii=where(logdata.footprint_id eq ids[i])
               if (ii[0] eq -1) then begin
                   Message, 'Could not find footprint id'
               endif
               idx_log[i]=ii
           endfor

           igbp=logdata(idx_log).igbp_index
       endif
   endif

   lat=dataG.latitude
   sza=dataG.solar_zenith_1

   xco2_tr=dataG.reported_xtarg_true
   xco2_rt=dataG.xtarg_ret
   xco2_sig=sqrt(dataG.xtarg_ret_var)

   ;; Scale true xco2 to ppm if needed
   if (xco2_tr[0] lt 1.0) then begin
       xco2_tr = xco2_tr * 1.e6
   endif
   
   ;; Scale retrieved xco2 to ppm if needed
   if (xco2_rt[0] lt 1.0) then begin
       xco2_rt = xco2_rt * 1.e6
       xco2_sig = xco2_sig * 1e6
   endif

   ; Get mean and find nearest value divisible by 5 to make the scale
   ; look better
   mean_xco2 = mean([xco2_tr, xco2_rt])
   for chk_down = fix(mean_xco2), fix(mean_xco2)-5, -1 do begin
       if (chk_down mod 5) eq 0 then $
         break
   endfor
   for chk_up = fix(mean_xco2), fix(mean_xco2)+5 do begin
       if (chk_up mod 5) eq 0 then $
         break
   endfor
   if abs(mean_xco2 - chk_down) lt abs(mean_xco2 - chk_up) then $
     mean_xco2 = chk_down $
   else $
     mean_xco2 = chk_up

   xco2_range = [mean_xco2 - 10, mean_xco2 + 10]
   
   aod_tr=dataG.reported_cont_oceanic_aod
   aod_tr_ice=dataG.reported_ice_aod
   aod_tr_h2o=dataG.reported_water_aod
   aod_tot_tr=aod_tr+aod_tr_ice+aod_tr_h2o
   
   aod_tot_rt=dataG.ret_all_total_aod
   aod_low_rt=dataG.ret_all_low_aod
   aod_mid_rt=dataG.ret_all_mid_aod
   aod_hi_rt=dataG.ret_all_high_aod
   aod_tot_err=aod_tot_rt-aod_tot_tr
   
   where_true_psurf = where(tag_names(dataG) eq "TRUE_PSURF")
   if where_true_psurf[0] ne -1 then $
     p_tr=dataG.true_psurf/100

   where_ap_psurf = where(tag_names(dataG) eq "AP_PSURF")
   if where_ap_psurf[0] ne -1 then $
     p_ap=dataG.ap_psurf/100

   where_true_psurf = where(tag_names(dataG) eq "RET_PSURF")
   if where_true_psurf[0] ne -1 then $
     p_rt=dataG.ret_psurf/100

   if n_elements(p_rt) gt 0 and n_elements(p_tr) gt 0 then $
     p_err=p_rt-p_tr

   if n_elements(p_rt) gt 0 and n_elements(p_ap) gt 0 then $
     p_rtap=p_rt-p_ap

   wndspd_tr=dataG.true_windspeed

   alb1_tr=dataG.true_albedo_1
   alb2_tr=dataG.true_albedo_2
   alb3_tr=dataG.true_albedo_3
   alb1_rt=dataG.ret_albedo_1
   alb2_rt=dataG.ret_albedo_2
   alb3_rt=dataG.ret_albedo_3
   alb1_err=alb1_rt-alb1_tr
   alb2_err=alb2_rt-alb2_tr
   alb3_err=alb3_rt-alb3_tr

   snr1=(statsG.beg_mean_snr_1+statsG.end_mean_snr_1)/2.
   snr2=(statsG.beg_mean_snr_2+statsG.end_mean_snr_2)/2.
   snr3=(statsG.beg_mean_snr_3+statsG.end_mean_snr_3)/2.
   chi2_1=statsG.chi2_1
   chi2_2=statsG.chi2_2
   chi2_3=statsG.chi2_3
   dsig=statsG.d_sigma_sq
   num_it=statsG.iterations
   rmsrel_1=statsG.rms_rel_1
   rmsrel_2=statsG.rms_rel_2
   rmsrel_3=statsG.rms_rel_3
   
   xco2_err=xco2_rt-xco2_tr
   tmp=moment(xco2_err)
   mn_err=tmp[0]
   std=sqrt(tmp[1])

   ;; Allow any further errors to be handled by default mechanisims
   Catch, /CANCEL

   ;; Try increasing threhsold for chi2 if no data is found in order
   ;; to ensure that something is plotted
   idx = -1
   while idx[0] eq -1 and max_chi2 lt 5.0 do begin
       ;; Combination of quality criteria from both Brian & Chris
       idx=where(chi2_1/chi2_2 lt max_rat12 and chi2_3/chi2_2 lt max_rat32 and $
                 chi2_1 lt max_chi2 and chi2_2 lt max_chi2 and chi2_3 lt max_chi2 and $
                 aod_tot_rt lt max_aod_rt and aod_tot_tr lt max_aod_tr and  $
                 xco2_sig lt max_xco2_sig, ngood)
       
       max_chi2  = max_chi2  + 1.0
       max_rat12 = max_rat12 + 1.0
       max_rat32 = max_rat32 + 1.0
   end

   if (ngood lt 2) then begin
       Message, /CONTINUE, 'Not enough data meets thresholds for files: ' + summary_file + ', ' + stats_file
       return
   endif

   print, 'Creating plot: ', out_filename

   ;; Values that will be reset at the end
   orig_p_font  = !P.font
   orig_p_thick = !P.thick
   orig_x_thick = !x.thick
   orig_y_thick = !y.thick
   orig_x_style = !x.style
   orig_y_style = !y.style

   orig_p_region   = !P.region
   orig_p_position = !P.position

   orig_p_color = !P.color
   orig_p_background = !P.background

   ;; Save the current color table.
   TVLCT, savedR, savedG, savedB, /GET

   ;; Initialize plots after sorting for valid data otherwise will
   ;; have empty plot files
   SetPS,/landscape,/color,hOffset=0.6,vOffset=0.6,/full
   ;InitColorBoss
   color=GetStandardPalette(/paper,/set,/darkblue,/blue,/darkred,/red,/lightgrey)
   device, filename=out_filename, /helvetica
   !P.thick=2
   !x.thick=2
   !y.thick=2
   !x.style=1
   !y.style=1

   qualityLabel='Precision < '+string(max_xco2_sig,'(f3.1)')+',  !Mc!X!U2!N (each band) < '+string(max_chi2,'(f3.1)')+'!c !Mc!X!U2!N (O!D2!N)/!Mc!X!U2!N (WCO!D2!N) < '+string(max_rat12,'(f3.1)')+', !Mc!X!U2!N (SCO!D2!N)/!Mc!X!U2!N (WCO!D2!N) < '+string(max_rat32,'(f3.1)')+'!c True Total AOD < '+string(max_aod_tr,'(f3.1)')+', Retrieved Total AOD < '+string(max_aod_rt,'(f3.1)')
  
   tmp=moment(xco2_err[idx])
   mn_gd=tmp[0]
   std_gd=sqrt(tmp[1])
   
   symbols,2,0.22
   pg=0

   ;; XCO2 (ret) vs XCO2 (truth), precision plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'Retrieved vs. True XCO!D2!N, Precision', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, xco2_tr,xco2_rt,psym=8,/noerase, $
         xrange=xco2_range,yrange=xco2_range, $
         xtitle='True XCO!D2!N / ppm',ytitle='Retrieved XCO!D2!N / ppm'
   oplot, xco2_tr[idx],xco2_rt[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[!y.cRange[0],!y.cRange[1]],linestyle=0
   oplot, [!x.cRange[0]+1.,!x.cRange[1]+1.],[!y.cRange[0]-1.,!y.cRange[1]-1.],linestyle=1
   oplot, [!x.cRange[0]-1.,!x.cRange[1]-1.],[!y.cRange[0]+1.,!y.cRange[1]+1.],linestyle=1

   ;xx=regress(xco2_tr[idx],xco2_rt[idx],correlation=cf,yfit=line)
   ;oplot,xco2_tr[idx],line

   Array_NextDivision
   plot, lat, xco2_tr,psym=8,/noerase,yrange=xco2_range,/nodata, $
         ytitle='Retrieved & True (red) XCO!D2!N / ppm', $
         xtitle='Latitude',xrange=[-90,90],subtitle='blue: '+qualityLabel
   symbols,2,0.25
   oplot, lat, xco2_tr,psym=8,color=color.red
   symbols,2,0.22
   oplot, lat, xco2_rt,psym=8
   oplot, lat[idx],xco2_rt[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, xco2_tr,xco2_err,psym=8,/noerase,xrange=xco2_range, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='True XCO!D2!N / ppm'
   oplot, xco2_tr[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Latitude',xrange=[-90,90], $
         subtitle='XCO!D2!N Mean Error = '+string(mn_gd,'(f5.2)')+'+/-'+string(std_gd,'(f4.2)')+'!c Number Good:  '+string(ngood,'(i5)')+' of '+string(ncases,'(i5)')
   oplot, lat[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, xco2_sig,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='XCO!D2!N Precision / ppm'
   oplot, xco2_sig[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,xco2_sig,psym=8,/noerase, $
         xtitle='Latitude',ytitle='XCO!D2!N Precision / ppm',xrange=[-90,90]
   oplot, lat[idx],xco2_sig[idx],psym=1,color=color.blue
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   Array_NextDivision

   ; sza, iterations plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'SZA, Number of Iterations', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, sza,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Solar Zenith Angle'
   oplot, sza[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,sza,psym=8,/noerase, $
         xtitle='Latitude',ytitle='Solar Zenith Angle',xrange=[-90,90]
   oplot, lat[idx],sza[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, num_it,xco2_err,psym=8,/noerase,xRange=[0,10],xTicks=5,xMinor=2, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Number of Iterations'
   oplot, num_it[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,num_it,psym=8,/noerase,yRange=[0,10], $
         xtitle='Latitude',ytitle='Number of Iterations',xrange=[-90,90]
   oplot, lat[idx],num_it,psym=1,color=color.blue
   Array_NextDivision
   Array_NextDivision
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   Array_NextDivision
   
   ; chi-squared plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'Individual Band !Mc!X!U2!N', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, chi2_1,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='!Mc!X!U2!N (A-Band)'
   oplot, chi2_1[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,chi2_1,psym=8,/noerase, $
         xtitle='Latitude',ytitle='!Mc!X!U2!N (A-Band)',xrange=[-90,90]
   oplot, lat[idx],chi2_1[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, chi2_2,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='!Mc!X!U2!N (Weak CO!D2!N)'
   oplot, chi2_2[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,chi2_2,psym=8,/noerase, $
         xtitle='Latitude',ytitle='!Mc!X!U2!N (Weak CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],chi2_2[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, chi2_3,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='!Mc!X!U2!N (Strong CO!D2!N)'
   oplot, chi2_3[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,chi2_3,psym=8,/noerase, $
         xtitle='Latitude',ytitle='!Mc!X!U2!N (Strong CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],chi2_3[idx],psym=1,color=color.blue
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   Array_NextDivision
   
   ; more chi-squared plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'Total !Mc!X!U2!N and Band Ratios', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, (chi2_1+chi2_2+chi2_3)/3.,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Total Measurement !Mc!X!U2!N'
   oplot, (chi2_1[idx]+chi2_2[idx]+chi2_3[idx])/3.,xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,(chi2_1+chi2_2+chi2_3)/3.,psym=8,/noerase, $
         xtitle='Latitude',ytitle='Total Measurement !Mc!X!U2!N',xrange=[-90,90]
   oplot, lat[idx],(chi2_1[idx]+chi2_2[idx]+chi2_3[idx])/3.,psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, chi2_1/chi2_2,xco2_err,psym=8,/noerase,xrange=[0.,2.], $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='!Mc!X!U2!N (A-Band) / !Mc!X!U2!N (Weak CO!D2!N)'
   oplot, chi2_1[idx]/chi2_2[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,chi2_1/chi2_2,psym=8,/noerase,yrange=[0.,2.], $
         xtitle='Latitude',ytitle='!Mc!X!U2!N (A-Band) / !Mc!X!U2!N (Weak CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],chi2_1[idx]/chi2_2[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, chi2_3/chi2_2,xco2_err,psym=8,/noerase,xrange=[0.,2.], $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='!Mc!X!U2!N (Strong CO!D2!N) / !Mc!X!U2!N (Weak CO!D2!N)'
   oplot, chi2_3[idx]/chi2_2[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,chi2_3/chi2_2,psym=8,/noerase,yrange=[0.,2.], $
         xtitle='Latitude',ytitle='!Mc!X!U2!N (Strong CO!D2!N) / !Mc!X!U2!N (Weak CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],chi2_3[idx]/chi2_2[idx],psym=1,color=color.blue
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   Array_NextDivision

   ; rms_rel plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'Individual Band Relative RMS', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, rmsrel_1,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Relative RMS (A-Band)'
   oplot, rmsrel_1[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,rmsrel_1,psym=8,/noerase, $
         xtitle='Latitude',ytitle='Relative RMS (A-Band)',xrange=[-90,90]
   oplot, lat[idx],rmsrel_1[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, rmsrel_2,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Relative RMS (Weak CO!D2!N)'
   oplot, rmsrel_2[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,rmsrel_2,psym=8,/noerase, $
         xtitle='Latitude',ytitle='Relative RMS (Weak CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],rmsrel_2[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, rmsrel_3,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Relative RMS (Strong CO!D2!N)'
   oplot, rmsrel_3[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,rmsrel_3,psym=8,/noerase, $
         xtitle='Latitude',ytitle='Relative RMS (Strong CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],rmsrel_3[idx],psym=1,color=color.blue
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   Array_NextDivision
   
   ; albedo plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'Individual Band Albedos', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   igd=where(alb1_rt gt -999)
   plot, alb1_rt[igd],xco2_err[igd],psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retrieved Albedo (A-Band)'
   oplot, alb1_rt[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat[igd],alb1_rt[igd],psym=8,/noerase, $
         xtitle='Latitude',ytitle='Retrieved Albedo (A-Band)',xrange=[-90,90]
   oplot, lat[idx],alb1_rt[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, alb2_rt[igd],xco2_err[igd],psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retrieved Albedo (Weak CO!D2!N)'
   oplot, alb2_rt[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat[igd],alb2_rt[igd],psym=8,/noerase, $
         xtitle='Latitude',ytitle='Retrieved Albedo (Weak CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],alb2_rt[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $      
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, alb3_rt[igd],xco2_err[igd],psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retrieved Albedo (Strong CO!D2!N)'
   oplot, alb3_rt[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat[igd],alb3_rt[igd],psym=8,/noerase, $
         xtitle='Latitude',ytitle='Retrieved Albedo (Strong CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],alb3_rt[idx],psym=1,color=color.blue
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   Array_NextDivision
   
   ; more albedo plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'Retrieved-True Albedos', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   igd=where(alb1_rt gt -999)
   plot, alb1_err[igd],xco2_err[igd],psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retvd-True Albedo (A-Band)'
   oplot, alb1_err[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   oplot, [0,0],[!y.cRange[0],!y.cRange[1]],linestyle=0
   Array_NextDivision
   plot, lat[igd],alb1_err[igd],psym=8,/noerase, $
         xtitle='Latitude',ytitle='Retvd-True Albedo (A-Band)',xrange=[-90,90]
   oplot, lat[idx],alb1_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[0,0], linestyle=0
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, alb2_err[igd],xco2_err[igd],psym=8,/noerase, $
         xRange=[-0.12,0.04],xTicks=4,xMinor=4, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retvd-True Albedo (Weak CO!D2!N)'
   oplot, alb2_err[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   oplot, [0,0],[!y.cRange[0],!y.cRange[1]],linestyle=0
   Array_NextDivision
   plot, lat[igd],alb2_err[igd],psym=8,/noerase, $
         xtitle='Latitude',ytitle='Retvd-True Albedo (Weak CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],alb2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[0,0], linestyle=0
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $      
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, alb3_err[igd],xco2_err[igd],psym=8,/noerase, $
         xRange=[-0.12,0.04],xTicks=4,xMinor=4, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retvd-True Albedo (Strong CO!D2!N)'
   oplot, alb3_err[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   oplot, [0,0],[!y.cRange[0],!y.cRange[1]],linestyle=0
   Array_NextDivision
   plot, lat[igd],alb3_err[igd],psym=8,/noerase, $
         xtitle='Latitude',ytitle='Retvd-True Albedo (Strong CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],alb3_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[0,0], linestyle=0
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   Array_NextDivision
   
   ; and even more albedo plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'Retrieved/True Albedos', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   igd=where(alb1_rt gt -999)
   plot, alb1_rt[igd]/alb1_tr[igd],xco2_err[igd],psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retvd/True Albedo (A-Band)'
   oplot, alb1_rt[idx]/alb1_tr[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   oplot, [1,1],[!y.cRange[0],!y.cRange[1]],linestyle=0
   Array_NextDivision
   plot, lat[igd],alb1_rt[igd]/alb1_tr[igd],psym=8,/noerase, $
         xtitle='Latitude',ytitle='Retvd/True Albedo (A-Band)',xrange=[-90,90]
   oplot, lat[idx],alb1_rt[idx]/alb1_tr[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1,1], linestyle=0
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, alb2_rt[igd]/alb2_tr[igd],xco2_err[igd],psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retvd/True Albedo (Weak CO!D2!N)'
   oplot, alb2_rt[idx]/alb2_tr[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   oplot, [1,1],[!y.cRange[0],!y.cRange[1]],linestyle=0
   Array_NextDivision
   plot, lat[igd],alb2_rt[igd]/alb2_tr[igd],psym=8,/noerase, $
         xtitle='Latitude',ytitle='Retvd/True Albedo (Weak CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],alb2_rt[idx]/alb2_tr[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1,1], linestyle=0
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $      
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, alb3_rt[igd]/alb3_tr[igd],xco2_err[igd],psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retvd/True Albedo (Strong CO!D2!N)'
   oplot, alb3_rt[idx]/alb3_tr[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   oplot, [1,1],[!y.cRange[0],!y.cRange[1]],linestyle=0
   Array_NextDivision
   plot, lat[igd],alb3_rt[igd]/alb3_tr[igd],psym=8,/noerase, $
         xtitle='Latitude',ytitle='Retvd/True Albedo (Strong CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],alb3_rt[idx]/alb3_tr[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1,1], linestyle=0
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   Array_NextDivision

   ; SNR plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'Individual Band SNR', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, snr1,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='SNR (A-Band)'
   oplot, snr1[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,snr1,psym=8,/noerase, $
         xtitle='Latitude',ytitle='SNR (A-Band)',xrange=[-90,90]
   oplot, lat[idx],snr1[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, snr2,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='SNR (Weak CO!D2!N)'
   oplot, snr2[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,snr2,psym=8,/noerase, $
         xtitle='Latitude',ytitle='SNR (Weak CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],snr2[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, snr3,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='SNR (Strong CO!D2!N)'
   oplot, snr3[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,snr3,psym=8,/noerase, $
         xtitle='Latitude',ytitle='SNR (Strong CO!D2!N)',xrange=[-90,90]
   oplot, lat[idx],snr3[idx],psym=1,color=color.blue
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm

   if n_elements(p_tr) gt 0 or n_elements(p_err) gt 0 or n_elements(p_rtap) gt 0 then begin
       Array_NextDivision
       
       ;; pressure plots
       Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                          title=plotTitle+'Surface Pressure', $
                          xMarginRatio=[0.08,0.05,0.05,0.02], $
                          yMarginRatio=[0.02,0.10]
   endif

   if n_elements(p_tr) gt 0 then begin
       Array_NewDivision, nx=1,ny=2,/byPosition, $
                          xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
       plot, p_tr,xco2_err,psym=8,/noerase, $
             ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='True P!DSurf!N / hPa'
       oplot, p_tr[idx],xco2_err[idx],psym=1,color=color.blue
       oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
       oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1

       Array_NextDivision
       plot, lat,p_tr,psym=8,/noerase, $
             xtitle='Latitude',ytitle='True P!DSurf!N / hPa',xrange=[-90,90]
       oplot, lat[idx],p_tr[idx],psym=1,color=color.blue
   endif

   if n_elements(p_err) gt 0 then begin
       Array_NextDivision
       Array_NewDivision, nx=1,ny=2,/byPosition, $
                          xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
       plot, p_err,xco2_err,psym=8,/noerase, $
             ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='P!DSurf!N Error (Rtvd-True) / hPa'
       oplot, p_err[idx],xco2_err[idx],psym=1,color=color.blue
       oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
       oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
       oplot, [0,0],[!y.cRange[0],!y.cRange[1]],linestyle=0
       Array_NextDivision
       plot, lat,p_err,psym=8,/noerase, $
             xtitle='Latitude',ytitle='P!DSurf!N Error (Rtvd-True) / hPa',xrange=[-90,90]
       oplot, lat[idx],p_err[idx],psym=1,color=color.blue
       oplot, [!x.cRange[0],!x.cRange[1]],[0,0],linestyle=0
   endif

   if n_elements(p_rtap) gt 0 then begin
       Array_NextDivision
       Array_NewDivision, nx=1,ny=2,/byPosition, $
                          xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
       plot, p_rtap,xco2_err,psym=8,/noerase, $
             ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='P!DSurf!N (Rtvd-apriori) / hPa'
       oplot, p_rtap[idx],xco2_err[idx],psym=1,color=color.blue
       oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
       oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
       oplot, [0,0],[!y.cRange[0],!y.cRange[1]],linestyle=0
       Array_NextDivision
       plot, lat,p_rtap,psym=8,/noerase, $
             xtitle='Latitude',ytitle='P!DSurf!N (Rtvd-apriori) / hPa',xrange=[-90,90]
       oplot, lat[idx],p_rtap[idx],psym=1,color=color.blue
       oplot, [!x.cRange[0],!x.cRange[1]],[0,0],linestyle=0
       pg=pg+1
   endif

   if n_elements(p_tr) gt 0 or n_elements(p_err) gt 0 or n_elements(p_rtap) gt 0 then begin
       xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   endif

   Array_NextDivision
   
   ; Aerosol Optical Depth plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'True Aerosol, Cloud Ice + H!D2!NO AOD', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, aod_tr,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='True Aerosol AOD'
   oplot, aod_tr[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,aod_tr,psym=8,/noerase, $
         xtitle='Latitude',ytitle='True Aerosol AOD',xrange=[-90,90]
   oplot, lat[idx],aod_tr[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, aod_tr_ice,xco2_err,psym=8,/noerase,xRange=[0.,0.3], $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='True Cloud Ice AOD'
   oplot, aod_tr_ice[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,aod_tr_ice,psym=8,/noerase,yRange=[0.,0.3], $
         xtitle='Latitude',ytitle='True Cloud Ice AOD',xrange=[-90,90]
   oplot, lat[idx],aod_tr_ice[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, aod_tr_h2o,xco2_err,psym=8,/noerase,xRange=[0.,1.0], $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='True Cloud Liquid H!D2!NO AOD'
   oplot, aod_tr_h2o[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,aod_tr_h2o,psym=8,/noerase,yRange=[0.,1.0], $
         xtitle='Latitude',ytitle='True Cloud Liquid H!D2!NO AOD',xrange=[-90,90]
   oplot, lat[idx],aod_tr_h2o[idx],psym=1,color=color.blue
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   Array_NextDivision

   ; more Aerosol Optical Depth plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'True Total AOD and Error', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, aod_tot_tr,xco2_err,psym=8,/noerase,xRange=[0.,1.5], $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='True Total AOD'
   oplot, aod_tot_tr[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,aod_tot_tr,psym=8,/noerase,yRange=[0.,1.5], $
         xtitle='Latitude',ytitle='True Total AOD',xrange=[-90,90]
   oplot, lat[idx],aod_tot_tr[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, aod_tot_err,xco2_err,psym=8,/noerase,xRange=[-1.0,1.0], $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Total AOD Error (Rtvd-True)'
   oplot, aod_tot_err[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   oplot, [0,0],[!y.cRange[0],!y.cRange[1]],linestyle=0
   Array_NextDivision
   plot, lat,aod_tot_err,psym=8,/noerase,yRange=[-1.0,1.0], $
         xtitle='Latitude',ytitle='Total AOD Error (Rtvd-True)',xrange=[-90,90]
   oplot, lat[idx],aod_tot_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[0,0],linestyle=0
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, aod_tot_rt,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retrieved Total AOD'
   oplot, aod_tot_rt[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, aod_tot_tr,aod_tot_rt,psym=8,/noerase,xRange=[0.,1.5],yRange=[0.,1.5], $
         ytitle='Retrieved Total AOD',xtitle='True Total AOD'
   oplot, aod_tot_tr[idx],aod_tot_rt[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[!y.cRange[0],!y.cRange[1]],linestyle=0
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   Array_NextDivision
   
   ; even more Aerosol Optical Depth plots
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'Low, Mid, High Retrieved AOD', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, aod_low_rt,xco2_err,psym=8,/noerase,xRange=[0.,0.8], $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retrieved Low AOD'
   oplot, aod_low_rt[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,aod_low_rt,psym=8,/noerase,yRange=[0.,0.8], $
         xtitle='Latitude',ytitle='Retrieved Low AOD',xrange=[-90,90]
   oplot, lat[idx],aod_low_rt[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, aod_mid_rt,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retrieved Mid AOD'
   oplot, aod_mid_rt[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,aod_mid_rt,psym=8,/noerase, $
         xtitle='Latitude',ytitle='Retrieved Mid AOD',xrange=[-90,90]
   oplot, lat[idx],aod_mid_rt[idx],psym=1,color=color.blue
   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, aod_hi_rt,xco2_err,psym=8,/noerase,xRange=[0.,0.3], $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='Retrieved High AOD'
   oplot, aod_hi_rt[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,aod_hi_rt,psym=8,/noerase,yRange=[0.,0.3], $
         xtitle='Latitude',ytitle='Retrieved High AOD',xrange=[-90,90]
   oplot, lat[idx],aod_hi_rt[idx],psym=1,color=color.blue
   Array_NextDivision
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   
   ; Land class, windspeed
   Array_NextDivision
   Array_NewDivision, nx=3,ny=1,/byPosition,/reset, $
                      title=plotTitle+'Surface Type, Windspeed', $
                      xMarginRatio=[0.08,0.05,0.05,0.02], $
                      yMarginRatio=[0.02,0.10]

   if n_elements(igbp) gt 0 then begin
       Array_NewDivision, nx=1,ny=2,/byPosition, $
                          xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
       plot, igbp,xco2_err,psym=8,/noerase,xRange=[0,18], $
             ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='IGBP Surface Type'
       oplot, igbp[idx],xco2_err[idx],psym=1,color=color.blue
       oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
       oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1

       Array_NextDivision
       plot, lat,igbp,psym=8,/noerase,yRange=[0,18], $
             xtitle='Latitude',ytitle='IGBP Surface Type',xrange=[-90,90]
       oplot, lat[idx],igbp[idx],psym=1,color=color.blue
   endif

   Array_NextDivision
   Array_NewDivision, nx=1,ny=2,/byPosition, $
                      xMarginRatio=[0.08,0.05],yMarginRatio=[0.02,0.12,0.04]
   plot, wndspd_tr,xco2_err,psym=8,/noerase, $
         ytitle='XCO!D2!N Error (Rtvd-True) / ppm',xtitle='True Windspeed'
   oplot, wndspd_tr[idx],xco2_err[idx],psym=1,color=color.blue
   oplot, [!x.cRange[0],!x.cRange[1]],[1.,1.],linestyle=1
   oplot, [!x.cRange[0],!x.cRange[1]],[-1.,-1.],linestyle=1
   Array_NextDivision
   plot, lat,wndspd_tr,psym=8,/noerase,yRange=[0,18], $
         xtitle='Latitude',ytitle='True Windspeed',xrange=[-90,90]
   oplot, lat[idx],wndspd_tr[idx],psym=1,color=color.blue
   Array_NextDivision
   pg=pg+1
   xyouts, 0.96,0.01,'Created '+strmid(systime(),8,2)+' '+strmid(systime(),4,3)+' '+strmid(systime(),20,4)+', page '+string(pg,'(i2)'),charsize=0.6,align=1.0,/norm
   Array_NextDivision

   device, /close_file

   ;; Reset global settings to values used on entry
   !P.font  = orig_p_font
   !P.thick = orig_p_thick
   !x.thick = orig_x_thick
   !y.thick = orig_y_thick
   !x.style = orig_x_style
   !y.style = orig_y_style

   !P.region = orig_p_region
   !P.position = orig_p_position

   !P.color = orig_p_color
   !P.background = orig_p_background

   tvlct, savedR, savedG, savedB

end

 
