PRO Reformat_Runlog_File, oldFile, newFile, BACKFORMAT=backFormat, CALTECH=caltech

   newRecord = { col1: '', $
                 runlab: '', $
                 iyr: 0L, $
                 iset:0L, $
                 zpdtim: 0.0D, $
                 oblat: 0.0D, $
                 oblon: 0.0D, $
                 zobs: 0.0D, $
                 asza: 0.0D, $
                 zenoff: 0.0D, $
                 opd: 0.0D, $
                 fovi: '', $
                 fovo: '', $
                 amal: '', $
                 ifirst: 0L, $
                 ilast: 0L, $
                 graw: 0.0D, $
                 possp: 0L, $
                 bytepw: 0L, $
                 zoff: 0.0D, $
                 snr: 0.0D, $
                 apf: '', $
                 tins: 0.0D, $
                 pins: 0.0D, $
                 hins: 0.0D, $
                 tout: 0.0D, $
                 pout: 0.0D, $
                 hout: 0.0D, $
                 lasf: 0.0D, $
                 wavtkr: 0.0D, $
                 sia: 0.0D, $
                 sis: 0.0D, $
                 aipl: 0.0D, $
                 spec_id: 0LL $
               }

   if keyword_set(caltech) then begin
       oldRecord = { col1: '', $
                     runlab: '', $
                     iyr: 0L, $
                     iset:0L, $
                     zpdtim: 0.0D, $
                     oblat: 0.0D, $
                     oblon: 0.0D, $
                     zobs: 0.0D, $
                     asza: 0.0D, $
                     zenoff: 0.0D, $
                     opd: 0.0D, $
                     fovi: '', $
                     fovo: '', $
                     amal: '', $
                     ifirst: 0L, $
                     ilast: 0L, $
                     graw: 0.0D, $
                     possp: 0L, $
                     bytepw: 0L, $
                     zoff: 0.0D, $
                     snr: 0.0D, $
                     apf: '', $
                     tins: 0.0D, $
                     pins: 0.0D, $
                     hins: 0.0D, $
                     tout: 0.0D, $
                     pout: 0.0D, $
                     hout: 0.0D, $
                     lasf: 0.0D, $
                     wavtkr: 0.0D, $
                     sia: 0.0D, $
                     sis: 0.0D, $
                     aipl: 0.0D $
       }
   endif else begin
       oldRecord = newRecord
   endelse

   if keyword_set(caltech) then $
     oldFormat = '(a1,a21,1x,2i4,f8.4,f8.3,f9.3,2f8.3,1x,f6.4,f7.2,3(1x,a5),2i8,1x,f14.11,i8,i3,1x,f5.3,i5,1x,a2,2(f6.1,f8.2,f5.1),f10.3,f7.0,2f6.1,f7.3)' $
   else $
     oldFormat = '(a1,a21,1x,2i4,f8.4,f8.3,f9.3,2f8.3,f7.0,f7.2,3(1x,a5),2i8,f15.11,i8,i3,1x,2f5.0,1x,a2,2(f6.0,f8.0,f5.0),f10.0,f7.0,2f6.1,f7.3,i18)'

   newFormat = '(a1,a35,1x,2i4,f8.4,f8.3,f9.3,2f8.3,1x,f6.4,f7.2,3(1x,a5),2i8,1x,f14.11,i8,i3,1x,f5.3,i5,1x,a2,2(f6.1,f8.2,f5.1),f10.3,f7.0,2f6.1,f7.3,i18)'

   IF N_Elements(newFile) LE 0 THEN $
     newFile = oldFile + ".new"

   OpenR, oldLun, oldFile, /GET_LUN, ERROR=openErr
   IF (openErr NE 0) THEN Message, !Error_State.msg, /NONAME

   OpenW, newLun, newFile, /GET_LUN, ERROR=openErr
   IF (openErr NE 0) THEN Message, !Error_State.msg, /NONAME

   headerStr=''
   ReadF, oldLun, headerStr

   PrintF, newLun, headerStr

   WHILE (eof(oldLun) NE 1) DO BEGIN
       IF Keyword_Set(backFormat) THEN BEGIN
           ReadF, oldLun, newRecord, FORMAT=newFormat 
           Print, 'New Format:'
           Print, newRecord, FORMAT=newFormat

           Copy_Struct_Contents, newRecord, oldRecord
       ENDIF ELSE BEGIN
           ReadF, oldLun, oldRecord, FORMAT=oldFormat
           Print, 'Old Format:'
           Print, oldRecord, FORMAT=oldFormat

           Copy_Struct_Contents, oldRecord, newRecord
       ENDELSE

       IF Keyword_Set(backFormat) THEN BEGIN
           us_pos  = strpos(oldRecord.runlab, '_')
           dot_pos = strpos(oldRecord.runlab, '.', /reverse_search)
           oldRecord.runlab = strmid(oldRecord.runlab, 0, us_pos) + strmid(oldRecord.runlab, dot_pos, 4) + ' '
       ENDIF ELSE BEGIN
           ;; Make strings that chop off leading so they will fit in allocated space
           newRecord.fovi = strmid(string(strcompress(newRecord.fovi, /REMOVE_ALL), format='(a6)'), 1)
           newRecord.fovo = strmid(string(strcompress(newRecord.fovo, /REMOVE_ALL), format='(a6)'), 1)
           newRecord.amal = strmid(string(strcompress(newRecord.amal, /REMOVE_ALL), format='(a6)'), 1)
       ENDELSE

       Print, ''

       IF Keyword_Set(backFormat) THEN BEGIN
           Print, 'Old Format:'
           PrintF, newLun, oldRecord, FORMAT=oldFormat
           Print, oldRecord, FORMAT=oldFormat
       ENDIF ELSE BEGIN
           Print, 'New Format:'
           PrintF, newLun, newRecord, FORMAT=newFormat
           Print, newRecord, FORMAT=newFormat
       ENDELSE
   ENDWHILE

   Free_Lun, oldLun
   Free_Lun, newLun

END

