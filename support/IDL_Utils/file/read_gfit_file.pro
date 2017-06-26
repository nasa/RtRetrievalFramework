FUNCTION Read_GFIT_File, gfitFilename

   recordChunkSize = 500

   OpenR, lun, gfitFilename, /GET_LUN, ERROR=openErr
   if (openErr ne 0) then Message, !Error_State.msg, /NONAME

   headerSkip = 0L
   numColumns = 0L
   ReadF, lun, headerSkip, numColumns

   ;; For some reason num_cols - 2 is what is reported. So fix here
   numColumns = numColumns + 2

   lineStr = ''
   for skipCount = 0, headerSkip-3 do begin
       ReadF, lun, lineStr
       ;print, 'header = ', lineStr
   endfor

   colNameStr = ''
   ReadF, lun, colNameStr
   columnNames = StrSplit(colNameStr, /EXTRACT)

   numColumns = min([numColumns, n_elements(columnNames)])

   columnNames = columnNames[0:numColumns-1]

   ;; Create matrix to store date
   gfitData = StrArr(numColumns, recordChunkSize)                   

   lineCount = 0
   WHILE (eof(lun) NE 1) DO BEGIN
       lineStr = ''
       ReadF, lun, lineStr
       lineData = StrSplit(lineStr, /EXTRACT)
       lineData = lineData[0:numColumns-1]

       numCurrRows = (size(gfitData))[2]
       if (lineCount gt numCurrRows-1) then begin
           ;print, 'Rebuilding to: ', numCurrRows + recordChunkSize
           gfitDataNew = StrArr(numColumns, numCurrRows + recordChunkSize)
           gfitDataNew[*, 0:lineCount-1] = gfitData[*, 0:lineCount-1]
           gfitData = gfitDataNew
       endif

       gfitData[*, lineCount] = lineData
       lineCount = lineCount + 1

   ENDWHILE

   Free_Lun, lun

   gfitData = gfitData[*, 0:lineCount-1]

   gfitObj = Obj_New('MATRIX_FILE')
   gfitObj->Set_All_Column_Labels, columnNames
   gfitObj->Set_Data, gfitData

   return, gfitObj
END
