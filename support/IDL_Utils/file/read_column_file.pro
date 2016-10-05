FUNCTION Read_Column_File, colFilename, HEADERSKIP=headerSkip, CHUNKSIZE=recordChunkSize, LIMIT=limit, VERBOSE=verbose

   IF NOT Keyword_Set(recordChunkSize) THEN $
     recordChunkSize = 500L

   OpenR, lun, colFilename, /GET_LUN, ERROR=openErr
   if (openErr ne 0) then Message, !Error_State.msg, /NONAME

   IF NOT Keyword_Set(headerSkip) THEN $
     headerSkip = 0L

   IF NOT Keyword_Set(limit) THEN $
     limit = -1L

   numColumns = 0L

   lineStr = ''
   for skipCount = 0, headerSkip-3 do begin
       ReadF, lun, lineStr
       ;print, 'header = ', lineStr
   endfor

   ;; Create matrix to store date
   ;;colData = StrArr(numColumns, recordChunkSize)                   

   lineCount = 0L
   WHILE ( (eof(lun) NE 1) AND (limit LE 0) OR (limit GT 0 AND lineCount LE limit) ) DO BEGIN
       lineStr = ''
       ReadF, lun, lineStr
       lineData = StrSplit(lineStr, /EXTRACT)

       numLineCols = N_Elements(lineData)
       newNumColumns = Max([numColumns, numLineCols])

       IF N_Elements(colData) GT 0 THEN $       
         numCurrRows = (size(colData))[2] $
       ELSE $
         numCurrRows = 0L

       IF (lineCount GT numCurrRows-1) OR (numColumns LT newNumColumns) THEN BEGIN
           if Keyword_Set(verbose) then $
             print, 'Rebuilding colData to size: ', newNumColumns, numCurrRows + recordChunkSize
           colDataNew = StrArr(newNumColumns, numCurrRows + recordChunkSize)

           IF N_Elements(colData) GT 0 THEN $
             colDataNew[0:numColumns-1, 0:lineCount-1] = colData[0:numColumns-1, 0:lineCount-1]

           colData = colDataNew
           numColumns = newNumColumns
       endif

       colData[0:numLineCols-1, lineCount] = lineData[0:numLineCols-1]
       lineCount = lineCount + 1

   ENDWHILE

   Free_Lun, lun

   colData = colData[*, 0:lineCount-1]

   colObj = Obj_New('MATRIX_FILE')
   colObj->Set_Data, colData

   return, colObj
END
