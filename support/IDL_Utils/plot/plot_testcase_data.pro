;; Loads a OCO L2 Type file containing testcase data where the first
;; column contains testcase names
PRO Plot_Testcase_Data, rmsInputFile, dataColumn, $
   ALL_FILTER=allFilter, $
   ALL_COLUMNS=allColumns, $
   SET_FILTER=setFilter, $
   SET_LINE_STYLES=streamsLineStyles, $
   SET_NAMES=setNames, $
   SET_SORT=setSort, $
   SET_COLORS=setColors, $
   AXIS_LABELS=axisLabels, $
   MIN_MAX_COLUMNS=minMaxCols, $
   SET_SHOWS_MIN_MAX=setShowsMinMax, $
   VARIANCE_COLUMN=varianceCol, $
   SET_SHOWS_VARIANCE=setShowsVariance, $
   _EXTRA=extra

   ;; Set up default values for keywords
   IF NOT Keyword_Set(setFilter) THEN $
     setFilter = [ '' ]

   IF N_Elements(setNames) LE 0 THEN $
     setNames = setFilter

   IF N_Elements(streamsLineStyles) LE 0 THEN $
     streamsLineStyles = IntArr( N_Elements(setFilter) ) ; All linestyle 0

   IF N_Elements(setColors) LE 0 THEN $
     setColors = IndGen( N_Elements(setNames) ) + 2

   IF N_Elements(allFilter) NE N_Elements(allColumns) THEN $
     Message, 'ALL_FILTER and ALL_COLUMNS keywords must be of same size'

   IF N_Elements(minMaxCols) GT 0 AND N_Elements(minMaxCols) NE 2 THEN $
     Message, 'MIN_MAX_COLUMNS must be a two index array [minIdx, maxIdx] if defined'

   IF N_Elements(setShowsMinMax) LE 0 THEN $
     setShowsMinMax = IntArr( N_Elements(setFilter ) ) + 1

   IF N_Elements(setShowsVariance) LE 0 THEN $
     setShowsVariance = IntArr( N_Elements(setFilter ) ) + 1

   ;; Load data from file
   tcDataObj = Obj_New('MATRIX_FILE', rmsInputFile)   
   tcDataMat = tcDataObj->Get_Data()

   ;; Filter all data for values matching the corresponding filter
   ;; strings and columns
   IF N_Elements(allFilter) GT 0 THEN BEGIN
       ;; Default to filter only on column 0
       IF N_Elements(allColumns) LE 0 THEN $
         allColumns = IntArr( N_Elements(allFilter) )

       FOR filtIdx = 0, N_Elements(allFilter) - 1 DO BEGIN
           filtCol = allColumns[filtIdx]
           filtRows = Where(StrMatch(tcDataMat[filtCol, *], '*' + allFilter[filtIdx] + '*') EQ 1)
           IF filtRows[0] NE -1 THEN $
             tcDataMat = tcDataMat[*, filtRows]
       ENDFOR
   ENDIF
   
   ;; Get run names column and specific data column
   runCol = tcDataMat[0, *]
   dataCol = Double( tcDataMat[dataColumn, *] )

   IF N_Elements(minMaxCols) GT 0 THEN BEGIN
       minCol = Double( tcDataMat[minMaxCols[0], *] )
       maxCol = Double( tcDataMat[minMaxCols[1], *] )
   ENDIF

   IF N_Elements(varianceCol) GT 0 THEN BEGIN
       varCol = Double( tcDataMat[varianceCol, *] )
   ENDIF
       
   IF NOT Keyword_Set(axisLabels) THEN $
     axisLabels = IndGen( N_Elements(dataCol) )
   
   Plot, axisLabels, axisLabels, YRANGE=[Min(dataCol), Max(dataCol)], COLOR=1, /NODATA, _EXTRA=extra
   
   FOR setIdx = 0, N_Elements(setFilter) - 1 DO BEGIN
       
       filtRows = Where(StrMatch(runCol, '*' + setFilter[setIdx] + '*') EQ 1)
       if filtRows[0] NE - 1 THEN BEGIN
           setData = dataCol[filtRows]
           setRuns = runCol[filtRows]

           IF N_Elements(minMaxCols) GT 0 THEN BEGIN
               setMin = minCol[filtRows]
               setMax = maxCol[filtRows]
           ENDIF

           IF N_Elements(varianceCol) GT 0 THEN BEGIN
               setVar = varCol[filtRows]
           ENDIF
       ENDIF ELSE BEGIN
           Continue
       ENDELSE
       
       IF N_Elements(setSort) GT 0 THEN BEGIN
           dataOrder = LonArr( N_Elements(setSort) )
           FOR sortIdx = 0, N_Elements(setSort) - 1 DO BEGIN
               labelLoc = Where(StrMatch(setRuns, '*' + setSort[sortIdx] + '*') EQ 1)
               dataOrder[sortIdx] = labelLoc[0]
           ENDFOR
           
           setData = setData[dataOrder]
           setRuns = setRuns[dataOrder]

           IF N_Elements(minMaxCols) GT 0 THEN BEGIN
               setMin = setMin[dataOrder]
               setMax = setMax[dataOrder]
           ENDIF

           IF N_Elements(varianceCol) GT 0 THEN BEGIN
               setVar = setVar[dataOrder]
           ENDIF
       ENDIF
       
       Print, setFilter[setIdx] + ': ', setRuns
       print, ''
       
       OPlot, axisLabels, setData, COLOR=setColors[setIdx], LINESTYLE=streamsLineStyles[setIdx], _EXTRA=extra

       IF N_Elements(minMaxCols) GT 0 AND setShowsMinMax[setIdx] THEN BEGIN
           Err_Plot, axisLabels, setMin, setMax, COLOR=setColors[setIdx]
       ENDIF

       IF N_Elements(varianceCol) GT 0 AND setShowsVariance[setIdx] THEN BEGIN
           stdDev = Sqrt( setVar )
           Err_Plot, axisLabels, setData-stdDev, setData+stdDev, COLOR=setColors[setIdx]
       ENDIF
       
   ENDFOR
   
   IF N_Elements(setNames) GT 0 THEN BEGIN
       txtColors = IntArr( N_Elements(setNames) ) + 1
       psymLeg = IntArr( N_Elements(setNames) ) + 0
       
       legend_dr, setNames, TEXTCOLORS=txtColors, COLORS=setColors, PSYM=psymLeg, $
                  CORNERS=corners, POS=[1, 1]
       
       xydims = [corners[2]-corners[0],corners[3]-corners[1]]
       chdim=[!d.x_ch_size/float(!d.x_size),!d.y_ch_size/float(!d.y_size)]
       pos = [!x.window[1]-chdim[0]-xydims[0], !y.window[1]-chdim[1]]
       
       
       legend_dr, setNames, TEXTCOLORS=txtColors, COLORS=setColors, PSYM=psymLeg, POS=pos
   ENDIF
   
END
