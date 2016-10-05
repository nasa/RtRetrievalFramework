;####################################################################
;
; Filename:        matrix_file__define.pro
;
; Class name:      Matrix_File
;
; Description:     Class representing OCO ASCII matrix data files
;
;########################### Change Log #############################
;
; Creator:              James McDuffie
; Creation date:        Sept. 18 2006
; Modifications:        See Subversion Notes
;
;####################################################################
; Copyright 2006, by the California Institute of Technology
; ALL RIGHTS RESERVED. United States Government Sponsorship
; acknowledged. Any commercial use must be negotiated with the Office
; of Technology Transfer at the California Institute of Technology.
;
; This software may be subject to U.S. export control laws and
; regulations. By accepting this document, the user agrees to comply
; with all applicable U.S. export laws and regulations. User has the
; responsibility to obtain export licenses, or other export authority
; as may be required before exporting such information to foreign
; countries or providing access to foreign persons.
;####################################################################


;####################################################################
;+
; Initializes the object with common header information
;-
;####################################################################

FUNCTION Matrix_File::Init, newFilename

   ;; Strings representing the beginning and end of the header
   self.beginHeaderStr = 'begin HEADER'
   self.endHeaderStr   = 'end HEADER'
   self.keywordSpacer  = '  '

   baseHeaderKeys = [ 'File_ID', 'File_Creation', 'File_Type' ]
   baseHeaderContents = StrArr(2, N_Elements(baseHeaderKeys))
   baseHeaderContents[0, *] = baseHeaderKeys

   self.p_header = Ptr_New(baseHeaderContents)
   self.p_data   = Ptr_New()
   self.p_formatFlag = Ptr_New()

   self->Set_Header_Keyword, 'File_ID', 'Unnamed Matrix Data'
   self->Set_Header_Keyword, 'File_Type', 'Matrix'

   IF N_Elements(newFilename) GT 0 THEN BEGIN
       self.filename = newFilename
       IF self.filename NE '' THEN $
         self->Read
   ENDIF ELSE BEGIN      
       self.filename = ''
   ENDELSE

   Return, 1B

END

;####################################################################
;+
; Clears up allocated memory
;-
;####################################################################

PRO Matrix_File::Cleanup

   Ptr_Free, self.p_header
   Ptr_Free, self.p_data

END

;####################################################################
;+
; Reads the file from disk
;-
;####################################################################

PRO Matrix_File::Read, newFilename

   IF N_Elements(newFilename) GT 0 THEN $
     IF newFilename NE '' THEN $
       self->Set_Filename, newFilename

   IF self.filename EQ '' THEN $
     Message, 'Filename not defined for class'

   ;; Error handling
   errStat = 0
   openErr = 0
;   Catch, errStat
   IF (errStat NE 0) THEN BEGIN
       ON_ERROR, 2
       Catch, /CANCEL
       IF (openErr EQ 0) THEN Free_Lun, lun
       Message, !Error_State.msg
   ENDIF

   ;; Open the file for read.
   OpenR, lun, self.filename, /GET_LUN, ERROR=openErr
   IF (openErr NE 0) THEN Message, !Error_State.msg, /NONAME

   ;; Parse header
   lineCount = 0
   finishedHeader = 0
   insideHeader = 0
   oldStyleFile = 0
   WHILE (NOT EOF(lun) AND NOT finishedHeader) DO BEGIN
       lineCount = lineCount + 1
       lineStr = ''
       ReadF, lun, lineStr

       ;; Trim off empty spaces from start and end of line
       lineStr = StrTrim(lineStr, 2)

       ;; Ignore empty lines and comments
       IF(lineStr EQ '') OR (StrMid(lineStr, 0, 1) EQ '#') THEN Continue

       ;; Check for Old style file.
       ;; The first line will have 3 dimensions of the file and the second
       ;; line will be the File_ID
       IF lineCount EQ 1 THEN BEGIN
           lineParts = StrSplit(lineStr, ' ', /EXTRACT)

           IF N_Elements(lineParts) EQ 3 THEN BEGIN
               oldStyleFile = 1

               self->Set_Header_Keyword, 'Num_Rows', lineParts[0]
               self->Set_Header_Keyword, 'Num_Columns', lineParts[1]
               skip = Fix( lineParts[2] )

               ;; Skip the specified number of lines. Assume first one
               ;; if it exists is the file ID
               IF skip GE 1 THEN BEGIN
                   ;; Read file ID line
                   fileID = ''
                   ReadF, lun, fileID
                   self->Set_Header_Keyword, 'File_ID', fileID

                   ;; Gobble up remaining skipped lines
                   FOR sCount = 1, skip-1 DO BEGIN
                       trashLine = ''
                       ReadF, lun, trashLine
                   ENDFOR

                   ;; Increment line counter to after skipped amount
                   lineCount = lineCount + skip
               ENDIF

               ;; Finished header since this was an old style file,
               ;; ready to start reading data
               finishedHeader = 1B
           ENDIF
       ENDIF

       ;; Found end of header
       IF(StrMid(lineStr, 0, StrLen(self.endHeaderStr)) EQ self.endHeaderStr) THEN BEGIN
           finishedHeader = 1
           insideHeader = 0
       ENDIF

       ;; Parsing lines inside the header
       IF insideHeader THEN BEGIN
           headerKeyVal = StrSplit(lineStr, '=', /EXTRACT, COUNT=kvCount)

           IF kvCount EQ 2 THEN BEGIN
               ;; Remove empty space
               keyword = StrTrim(headerKeyVal[0], 2)
               value   = StrTrim(headerKeyVal[1], 2)
          
               self->Set_Header_Keyword, keyword, value
           ENDIF
       ENDIF

       ;; Found start of header
       IF(StrMid(lineStr, 0, StrLen(self.beginHeaderStr)) EQ self.beginHeaderStr) THEN $
         insideHeader = 1

   ENDWHILE

   ;; Intialize the data structure with the newly read num rows and
   ;; columns header keywords
   if self->Header_Has_Keyword('Num_Rows') then begin
       numRows = self->Get_Header_Keyword('Num_Rows')
   endif else begin
       numRows = 0
   endelse
   
   if self->Header_Has_Keyword('Num_Columns') then begin
       numCols = self->Get_Header_Keyword('Num_Columns') 
   endif else if numRows gt 0 then begin
       numCols = 1 
   endif else begin
       numCols = 0
   endelse

   ;; Only try to read data if any exists, otherwise file just
   ;; consists of a header
   if numRows GT 0L and numCols GT 0L THEN BEGIN

       ;; Load values into new data structure
       IF oldStyleFile THEN BEGIN
           newData = DblArr( numCols, numRows )

           ;; No comments in old style files and it is not guaranteed that
           ;; each line will contain all parts of a row
           ReadF, lun, newData
       ENDIF ELSE BEGIN
           newData = StrArr( numCols, numRows )

           ;; Loop over each line since a line could possibly be a comment
           rowIndex = 0L
           columnsLeft = numCols
           columnStart = 0L
           WHILE (NOT EOF(lun)) DO BEGIN
               lineCount = lineCount + 1
               lineStr = ''
               ReadF, lun, lineStr

               ;; Trim off empty spaces from start and end of line
               lineStr = StrTrim(lineStr, 2)
               
               ;; Ignore empty lines and comments
               IF(lineStr EQ '') OR (StrMid(lineStr, 0, 1) EQ '#') THEN Continue

               newRow = StrSplit(lineStr, /EXTRACT)

               IF numRows LE rowIndex THEN $
                 break

               newData[columnStart:columnStart+N_Elements(newRow)-1, rowIndex] = newRow
               columnsLeft = columnsLeft - N_Elements(newRow)

               if columnsLeft eq 0L then begin
                   columnStart = 0L
                   columnsLeft = numCols
                   rowIndex = rowIndex + 1L
               endif else begin
                   columnStart = columnStart + N_Elements(newRow)
               endelse
           ENDWHILE

           IF rowIndex NE numRows THEN $
             Message, /CONTINUE, String( 'Number of rows read: ', rowIndex, $
                                         ' does not match expected number: ', $
                                         numRows, ' in file: ', self.filename, $
                                         FORMAT='(A0, I0, A0, I0, A0, A0)' )

           
           newData = newData[0L:numCols-1L, 0L:numRows-1L]
       ENDELSE

       self->Set_Data, newData
   ENDIF

   ;; Close the file
   Free_Lun, lun

END

;####################################################################
;+
; Writes the contents of the matrix data to a file
;-
;####################################################################

PRO Matrix_File::Write, newFilename, COLUMN_FORMAT=matrixFormatCodes

   IF N_Elements(newFilename) GT 0 THEN $
     IF newFilename NE '' THEN $
       self->Set_Filename, newFilename

   IF N_Elements(matrixFormatCodes) GT 0 THEN $
     IF N_Elements(matrixFormatCodes) NE self->Get_Num_Columns() THEN $
     Message, 'Format codes must be defined for all columns'

   ;; Error handling
   errStat = 0
   ;;Catch, errStat
   IF (errStat NE 0) THEN BEGIN
       ON_ERROR, 2
       Catch, /CANCEL
       IF (openErr EQ 0) THEN Free_Lun, lun
       Message, !Error_State.msg
   ENDIF

   ;; Set the file creation time if the keyword does not already exist
   ;; or is empty
   IF self->Header_Has_Keyword('File_Creation') THEN BEGIN
       IF self->Get_Header_Keyword('File_Creation') EQ '' THEN $
         self->Set_Header_Keyword, 'File_Creation', SysTime()
   ENDIF ELSE BEGIN
       self->Set_Header_Keyword, 'File_Creation', SysTime()
   ENDELSE

   ;; Open the file for write, and get a lun.
   OpenW, lun, newFilename, /GET_LUN, ERROR=openErr
   IF (openErr NE 0) THEN Message, !Error_State.msg, /NONAME
   
   ;; Begin header string
   PrintF, lun, self.beginHeaderStr

   ;; Write the file header keywords
   headerNames = self->Get_Header_Keyword_Names()
   FOR keyIdx = 0, N_Elements(headerNames)-1 DO BEGIN
       keyword = headerNames[keyIdx]
       value   = self->Get_Header_Keyword(keyword)
       
       PrintF, lun, self.keywordSpacer + keyword + ' = ' + value
   ENDFOR

   ;; End header string
   PrintF, lun, self.endHeaderStr

   ;; Retrieve matrix data size
   numRows = self->Get_Num_Rows()
   numCols = self->Get_Num_Columns()

   IF N_Elements(matrixFormatCodes) EQ 0 THEN BEGIN
       ;; Create column size formating strings
       labelFormatCodes = StrArr( numCols )
       matrixFormatCodes = StrArr( numCols )
       FOR colIdx = 0, numCols-1 DO BEGIN
           colPrec = self->Get_Column_Precision(colIdx, PRECTYPE=precType)
           precStr = String(colPrec, FORMAT='(I0)')
           
           labelWidth = 0
           IF self->Header_Has_Keyword('Labels') THEN $
             labelWidth = StrLen( self->Get_Column_Label(colIdx) )

           colWidth = colPrec
           IF precType EQ 'E' THEN $
             colWidth = colWidth + 8 $
           ELSE $
             colWidth = colWidth + 1

           colWidth = Max([colWidth, labelWidth + 1])
           widthStr = String(colWidth, FORMAT='(I0)')
           
           IF precType EQ 'E' THEN $
             matrixFormatCodes[colIdx] = precType + widthStr + '.' + precStr $
           ELSE $
             matrixFormatCodes[colIdx] = precType + widthStr

           labelFormatCodes[colIdx] = 'A' + widthStr

       ENDFOR
   ENDIF ELSE BEGIN
       labelFormatCodes = StrArr( numCols )
       FOR colIdx = 0, numCols-1 DO BEGIN
           labelFormatCodes[colIdx] = String(Strlen(String(0, format='('+matrixFormatCodes[colIdx]+')')), format='("A", I0)')
       ENDFOR
   ENDELSE

   labelFormat = '(' + StrJoin(labelFormatCodes, ',') + ')'
   matrixFormat = '(' + StrJoin(matrixFormatCodes, ',') + ')'

   ;; Write labels as comment before
   IF (self->Header_Has_Keyword('Labels')) THEN BEGIN
       colLabels = self->Get_All_Keyword_Items('Labels')
       PrintF, lun, '#' + String(colLabels, FORMAT=labelFormat)
   ENDIF

   ;; Write units as comments
   IF (self->Header_Has_Keyword('Units')) THEN BEGIN
       colUnits = '(' + self->Get_All_Keyword_Items('Units') + ')'
       PrintF, lun, '#' + String(colUnits, FORMAT=labelFormat)
   ENDIF

   ;; Write matrix contents
   FOR rowIdx = 0, numRows-1 DO BEGIN
       PrintF, lun, ' ' + String((*self.p_data)[*, rowIdx], FORMAT=matrixFormat)
   ENDFOR

   ;; Close file
   Free_Lun, lun

END

;####################################################################
;+
; Gets the filename of the data represented by the class
;-
;####################################################################

FUNCTION Matrix_File::Get_Filename
   Return, self.filename
END

;####################################################################
;+
; Sets the filename for the data represented by the class
;-
;####################################################################

PRO Matrix_File::Set_Filename, newFilename
   self.filename = newFilename
END

;####################################################################
;+
; Returns a header keyword's value stored it in the class
;-
;####################################################################

FUNCTION Matrix_File::Get_Header_Keyword, keyword

   ;; Return the stored value or an error if the value does not exist
   IF self->Header_Has_Keyword(keyword, INDEX=keywordIndex) THEN BEGIN
       Return, (*self.p_header)[ 1, keywordIndex ]
   ENDIF ELSE BEGIN
       Message, 'Header keyword: ' + keyword + ' not found in file: ' + $
                self.filename
   ENDELSE

END

;####################################################################
;+
; Takes a header keyword and value stores it in the class
;-
;####################################################################

PRO Matrix_File::Set_Header_Keyword, keyword, value

   ;; Create a new location for the keyword since it does not
   ;; already exist
   IF self->Header_Has_Keyword(keyword, INDEX=keywordIndex) THEN BEGIN
       (*self.p_header)[ 0, keywordIndex ] = keyword
       (*self.p_header)[ 1, keywordIndex ] = value      
   ENDIF ELSE BEGIN
       numKeywords = N_Elements( (*self.p_header)[0, *] )

       newHeader = StrArr(2, numKeywords + 1)

       newHeader[*, 0:numKeywords-1] = *self.p_header
       Ptr_Free, self.p_header

       newHeader[0, numKeywords] = keyword
       newHeader[1, numKeywords] = value

       self.p_Header = Ptr_New(newHeader)
   ENDELSE
   
END

;####################################################################
;+
; Returns the name of all header keywords
;-
;####################################################################

FUNCTION Matrix_File::Get_Header_Keyword_Names

   Return, Transpose( (*self.p_header)[0, *] )

END

;####################################################################
;+
; Returns whether or not the header has a particular keyword
;-
;####################################################################

FUNCTION Matrix_File::Header_Has_Keyword, keywordQuery, INDEX=keyIndex

   headerNames = self->Get_Header_Keyword_Names()

   whereKeyword = Where(StrUpCase(headerNames) EQ StrUpCase(keywordQuery))

   hasKeyword = whereKeyword[0] NE -1
   keyIndex = whereKeyword[0]

   Return, hasKeyword

END

;####################################################################
;+
; Returns a boolean for whether or not the class contains valid data
;-
;####################################################################

FUNCTION Matrix_File::Has_Valid_Data

   Return, Ptr_Valid(self.p_data)

END

;####################################################################
;+
; Returns all the data values stored by the class
;-
;####################################################################

FUNCTION Matrix_File::Get_Data, NOCASTDOUBLE=nocastdouble

   IF Ptr_Valid(self.p_data) THEN BEGIN

       if Keyword_Set(nocastdouble) then $
         castDouble = 0B $
       else $
         castDouble = 1B

       Catch, errStat
       IF (errStat NE 0) THEN BEGIN
           Catch, /CANCEL
           type_error:
           castDouble = 0B
       ENDIF

       On_IOError, type_error

       IF castDouble THEN BEGIN
           Return, Double( *self.p_data )
       ENDIF ELSE BEGIN
           Return, *self.p_data
       ENDELSE

   ENDIF ELSE BEGIN
       Message, 'Data not yet initialized by class.'
   ENDELSE

END

;####################################################################
;+
; Gets a column of data
;-
;####################################################################

FUNCTION Matrix_File::Get_Column_Data, columnId, NOCASTDOUBLE=nocastdouble

   colIndex = self->Get_Index_Of_Label(columnId)
   
   IF colIndex LT 0 OR colIndex GE self->Get_Num_Columns() THEN $
     Message, $
       String( 'Column index: ', colIndex, $
               ' out of range of number of columns: ', $
               self->Get_Num_Columns(), ' in file: ', self.filename, $
               FORMAT='(A0, I0, A0, I0, A0, A0)' )

   IF Ptr_Valid(self.p_data) THEN BEGIN

       if Keyword_Set(nocastdouble) then $
         castDouble = 0B $
       else $
         castDouble = 1B

       Catch, errStat
       IF (errStat NE 0) THEN BEGIN
           Catch, /CANCEL
           type_error:
           castDouble = 0B
       ENDIF

       On_IOError, type_error

       IF castDouble THEN $
         Return, Double( ((*self.p_data)[colIndex, *])[*] ) $
       ELSE $
         Return, ((*self.p_data)[colIndex, *])[*]

   ENDIF ELSE BEGIN
       Message, 'Data not yet initialized by class'
   ENDELSE

END

;####################################################################
;+
; Returns all the data values stored by the class in a struct array
;-
;####################################################################

FUNCTION Matrix_File::Get_Struct_Array

   IF Ptr_Valid(self.p_data) THEN BEGIN

       columnLabels = self->Get_All_Column_Labels()

       createCmd = 'elt_ = {'
       for colIndex = 0, self->Get_Num_Columns()-1 do begin
           colName = StrUpCase(StrTrim(columnLabels[colIndex]))

           if IS_Reserved_Word(colName) then begin
               colName = colName + '_COLUMN'
           endif
           
           precision = self->Get_Column_Precision(colIndex, PRECTYPE=precType)
           
           if precType eq 'I' then $
             initVal = '0L' $
           else if precType eq 'E' then $
             initVal = '0.0D0' $
           else $
             initVal = '" "'

           createCmd = createCmd + colName + ':' + initVal
           if (colIndex ne (self->Get_Num_Columns()-1)) then createCmd = createCmd + ', '
       endfor
       
       createCmd = createCmd + '}'
       errorResult = execute(createCmd)
       
       structArray = replicate(elt_, self->Get_Num_Rows())

       for rowIdx = 0, self->Get_Num_Rows()-1 do begin
           for colIdx = 0, self->Get_Num_Columns()-1 do begin
               structArray[rowIdx].(colIdx) = (*self.p_data)[colIdx, rowIdx]
           endfor
       endfor

       Return, structArray
   ENDIF ELSE BEGIN
       Message, 'Data not yet initialized by class.'
   ENDELSE

END



;####################################################################
;+
; Sets a new matrix as the class data
;-
;####################################################################

PRO Matrix_File::Set_Data, newData

   Ptr_Free, self.p_data

   self.p_data = Ptr_New( String(newData) )

   self->Set_Header_Keyword, 'Num_Rows', $
     String( self->Get_Num_Rows(), FORMAT='(I0)' )
   self->Set_Header_Keyword, 'Num_Columns', $
     String( self->Get_Num_Columns(), FORMAT='(I0)' )

END

;####################################################################
;+
; Returns the number of rows stored by the class
;-
;####################################################################

FUNCTION Matrix_File::Get_Num_Rows

   IF Ptr_Valid(self.p_data) THEN $
     Return, N_Elements( (*self.p_data)[0, *] ) $
   ELSE $
     Return, 0

END

;####################################################################
;+
; Returns the number of columns stored by the class
;-
;####################################################################

FUNCTION Matrix_File::Get_Num_Columns

   IF Ptr_Valid(self.p_data) THEN $
     Return, N_Elements( (*self.p_data)[*, 0] ) $
   ELSE $
     Return, 0

END

;####################################################################
;+
; Sets flag indicating which columns should have their precision
; tested. If 0 then the column is treated as a string.
;-
;####################################################################

PRO Matrix_File::Set_Column_Format_Flag, newFlags

   Ptr_Free, self.p_formatFlag

   self.p_formatFlag = Ptr_New( Fix(newFlags) )

END


;####################################################################
;+
; Returns the least amount of decimals needed in exponential notation
; to preserve the precision of the column
;-
;####################################################################

FUNCTION Matrix_File::Get_Column_Precision, colIndex, PRECTYPE=precType

   IF colIndex LT 0 OR colIndex GE self->Get_Num_Columns() THEN $
     Message, $
       String( 'Column index: ', colIndex, $
               ' out of range of number of columns: ', $
               self->Get_Num_Columns(), ' in file: ', self.filename, $
               FORMAT='(A0, I0, A0, I0, A0, A0)' )

   colData = (*self.p_data)[colIndex, *]

   IF Ptr_Valid(self.p_formatFlag) THEN BEGIN
       IF colIndex LT N_Elements(self.p_formatFlag) THEN BEGIN
           IF (*self.p_formatFlag)[colIndex] NE 1 THEN BEGIN
               bestPrecision = $
                 Max(Strlen(String(colData, FORMAT='(A0)')))
               precType = 'A'
               Return, bestPrecision
           ENDIF
       ENDIF
   ENDIF

   eps = (MachAr(/DOUBLE)).eps

   valid_num_loc = IntArr( N_Elements(colData) )
   FOR data_idx = 0, N_Elements(colData)-1 DO BEGIN
       valid_num_loc[data_idx] = Valid_Num( colData[data_idx] )
   ENDFOR

   where_num = Where(valid_num_loc eq 1, valid_count)

   IF valid_count eq N_Elements(colData) THEN BEGIN
       holdsPrecision = 1B
       bestPrecision = 20L
       currPrecision = bestPrecision
       WHILE (holdsPrecision AND currPrecision GE 0) DO BEGIN
           ;; Turn precision integer into string for formatting
           precStr = String(currPrecision, FORMAT='(I0)')

           ;; Convert column to string w/ current precision then back to
           ;; a double and compare against the original column data and
           ;; machine precision
           colCompare = Double(String(colData, FORMAT='(E0.'+precStr+')'))
           colPrecMeasure = $
             Abs( colCompare - colData  ) LT eps
           
           ;; Make sure all are less than machine noise in differents and
           ;; mark the current precision as best seen
           whereOne = Where(colPrecMeasure EQ 1, oneCount)
           IF oneCount EQ N_Elements(colData) THEN BEGIN
               bestPrecision = currPrecision
           ENDIF ELSE BEGIN
               holdsPrecision = 0B
           ENDELSE

           currPrecision = currPrecision - 1

           ; Make sure that item is really an int by comparing the conversion
           ; of the column data to ints with the real numbers
           intColMeasure = Abs( fix(colData) - colCompare  ) LT eps
           
       ENDWHILE

       whereTrueInt = Where(intColMeasure EQ 1, trueIntCount)

       IF trueIntCount EQ N_Elements(colData) THEN BEGIN
           bestPrecision = $
             Max(Strlen(String(colData, FORMAT='(I0)')))
           precType = 'I'
       ENDIF ELSE BEGIN
           IF bestPrecision EQ 0 THEN $
             bestPrecision = 2
           precType = 'E'
       ENDELSE
   ENDIF ELSE BEGIN ;; IF valid number
       bestPrecision = $
         Max(Strlen(String(colData, FORMAT='(A0)')))
       precType = 'A'
   ENDELSE

   Return, bestPrecision

END

;####################################################################
;+
; Returns all column label names
;-
;####################################################################

FUNCTION Matrix_File::Get_All_Column_Labels

   Return, self->Get_All_Keyword_Items('Labels')
     
END

;####################################################################
;+
; Sets all column label names
;-
;####################################################################

PRO Matrix_File::Set_All_Column_Labels, newLabels

   IF self->Get_Num_Columns() EQ 0 THEN BEGIN
       newLabels[*] = '"' + newLabels[*] + '"'
       labelsStr = StrJoin(newLabels, ' ')

       self->Set_Header_Keyword, 'Labels', labelsStr
   ENDIF ELSE BEGIN
       self->Set_All_Keyword_Items, 'Labels', newLabels
   ENDELSE
     
END


;####################################################################
;+
; Returns the label for a particular column index
;-
;####################################################################

FUNCTION Matrix_File::Get_Column_Label, colIndex  

   IF colIndex LT 0 OR colIndex GE self->Get_Num_Columns() THEN $
     Message, $
       String( 'Column index: ', colIndex, $
               ' out of range of number of columns: ', $
               self->Get_Num_Columns(), ' in file: ', self.filename, $
               FORMAT='(A0, I0, A0, I0, A0, A0)' )

   colLabels = self->Get_All_Keyword_Items('Labels')

   Return, colLabels[colIndex]
     
END

;####################################################################
;+
; Sets the label for a particular column index
;-
;####################################################################

PRO Matrix_File::Set_Column_Label, colIndex, newLabel
   
   IF colIndex LT 0 OR colIndex GE self->Get_Num_Columns() THEN $
     Message, $
       String( 'Column index: ', colIndex, $
               ' out of range of number of columns: ', $
               self->Get_Num_Columns(), ' in file: ', self.filename, $
               FORMAT='(A0, I0, A0, I0, A0, A0)' )

   colLabels = self->Get_All_Keyword_Items('Labels')

   colLabels[colIndex] = newLabel

   self->Set_All_Keyword_Items, 'Labels', colLabels

END

;####################################################################
;+
; Gets an index of a column name
;-
;####################################################################

FUNCTION Matrix_File::Get_Index_Of_Label, columnId

   colIndex = -1L
   IF Valid_Num(columnId) THEN BEGIN
       colIndex = columnId
   ENDIF ELSE BEGIN
       whereId = Where(StrLowCase(columnId) EQ StrLowCase(self->Get_All_Keyword_Items('Labels')))
       IF whereId[0] NE -1 THEN $
         colIndex = whereId[0] $
       ELSE $
         Message, 'Column name: ' + columnId + ' not found.'
   ENDELSE

   Return, colIndex
END

;####################################################################
;+
; Returns the unit for a particular column index
;-
;####################################################################

FUNCTION Matrix_File::Get_Column_Unit, colIndex  

   IF colIndex LT 0 OR colIndex GE self->Get_Num_Columns() THEN $
     Message, $
       String( 'Column index: ', colIndex, $
               ' out of range of number of columns: ', $
               self->Get_Num_Columns(), ' in file: ', self.filename, $
               FORMAT='(A0, I0, A0, I0, A0, A0)' )

   colUnits = self->Get_All_Keyword_Items('Units')

   Return, colUnits[colIndex]
     
END

;####################################################################
;+
; Sets the unit for a particular column index
;-
;####################################################################

PRO Matrix_File::Set_Column_Unit, colIndex, newUnit
   
   IF colIndex LT 0 OR colIndex GE self->Get_Num_Columns() THEN $
     Message, $
       String( 'Column index: ', colIndex, $
               ' out of range of number of columns: ', $
               self->Get_Num_Columns(), ' in file: ', self.filename, $
               FORMAT='(A0, I0, A0, I0, A0, A0)' )

   colUnits = self->Get_All_Keyword_Items('Units')

   colUnits[colIndex] = newUnit

   self->Set_All_Keyword_Items, 'Units', colUnits

END

;####################################################################
;+
; Retrieves all column labels
;-
;####################################################################

FUNCTION Matrix_File::Get_All_Keyword_Items, keyword

   colLabels = Replicate( 'NO_LABEL', self->Get_Num_Columns() )

   IF self->Header_Has_Keyword(keyword) THEN BEGIN
       labelsStr = self->Get_Header_Keyword(keyword)
       labelsArr = StrSplit(labelsStr, '(^"|"[ ]*"|"$)', /EXTRACT, /REGEX)
       numSetLabels = N_Elements(labelsArr)

       IF numSetLabels GT self->Get_Num_Columns() THEN $
         Message, 'Error extracting item list from "' + keyword + '" ' + $
                  'header keyword in file: ' + self.filename

       IF numSetLabels GT 0 THEN $
         colLabels[0:numSetLabels-1] = labelsArr
   ENDIF

   return, colLabels
   
END

;####################################################################
;+
; Sets all column labels
;-
;####################################################################

PRO Matrix_File::Set_All_Keyword_Items, keyword, newColumnLabels

   numExistingCols = self->Get_Num_Columns()
   numNewCols = N_Elements(newColumnLabels)

   IF numNewCols NE numExistingCols THEN $
     Message, /CONTINUE, $
       String( 'WARNING: Number of new column labels: ', numNewCols, $
               ' does not match number of data columns: ', $
               self->Get_Num_Columns(), $
               '. Using only those with indexes matching current columns. ', $
               ' in file: ', self.filename, $
               FORMAT='(A0, I0, A0, I0, A0, A0, A0)' )

   newColumnLabels[0:numExistingCols-1] = $
     '"' + newColumnLabels[0:numExistingCols-1] + '"'
   labelsStr = StrJoin(newColumnLabels[0:numExistingCols-1], ' ')

   self->Set_Header_Keyword, keyword, labelsStr

END

;####################################################################
;+
; Returns the header File_Creation keyword
;-
;####################################################################

FUNCTION Matrix_File::Get_File_Creation

   creation = self->Get_Header_Keyword('File_Creation')
   creationStr = StrSplit( creation, '"', ESCAPE='\', /EXTRACT )

   Return, creationStr

END

;####################################################################
;+
; Sets the header File_Creation keyword
;-
;####################################################################

PRO Matrix_File::Set_File_Creation, newCreation

   self->Set_Header_Keyword, 'File_Creation', '"' + newCreation + '"'

END

;####################################################################
;+
; Gets the header File_ID keyword
;-
;####################################################################

FUNCTION Matrix_File::Get_File_ID

   fileId = self->Get_Header_Keyword('File_ID')
   fileIdStr = StrSplit( fileId, '"', ESCAPE='\', /EXTRACT )

   Return, fileIdStr[0]

END

;####################################################################
;+
; Sets the header File_ID keyword
;-
;####################################################################

PRO Matrix_File::Set_File_ID, newFileID

   self->Set_Header_Keyword, 'File_ID', '"' + newFileID + '"'
END

;####################################################################
;+
; Gets the header File_Type keyword
;-
;####################################################################

FUNCTION Matrix_File::Get_File_Type

   fileType = self->Get_Header_Keyword('File_Type')
   fileTypeStr = StrSplit( fileType, '"', ESCAPE='\', /EXTRACT )

   Return, fileTypeStr[0]

END

;####################################################################
;+
; Sets the header File_Type keyword
;-
;####################################################################

PRO Matrix_File::Set_File_Type, newFileType

   self->Set_Header_Keyword, 'File_Type', '"' + newFileType + '"'

END


;####################################################################
;+
; Create structure defining attributes of the class
;-
;####################################################################

PRO Matrix_File__Define

   temp = { MATRIX_FILE,                 $
            filename        : '',        $
            p_header       : Ptr_New(),  $
            p_data         : Ptr_New(),  $
            p_formatFlag   : Ptr_New(),  $
            beginHeaderStr : '',         $
            endHeaderStr   : '',         $
            keywordSpacer  : ''          $
          }

END
