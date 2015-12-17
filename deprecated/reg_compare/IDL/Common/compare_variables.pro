;+
;;##################### TES Level 1B IDL Code ########################
;;
;;                      Tropospheric Emission
;;                          Spectrometer
;;                             (TES)
;;
;;                       Level 1B Subsystem
;;
;;                         Copyright 2004,
;;              California Institute of Technology.
;;                       ALL RIGHTS RESERVED.
;;           U.S. Government Sponsorship acknowledged.
;;
;;
;; ROUTINE NAME
;;    Compare_Variables
;;
;; FILE NAME
;;    compare_variables.pro
;;
;; ABSTRACT
;;    Compares two IDL variables
;; 
;; AUTHOR
;;    James McDuffie
;;
;; DESCRIPTION
;;    Compares two variables of the same type and returns the
;;    number of differences. An optional string called report can be
;;    used for getting difference information.

FUNCTION Compare_Variables, baseVal, $
                            compVal, $
                            name, $
                            checkAll, $
                            REPORT=report, $
                            _EXTRA=extra

   ;; Reset report
   report = ''

   ;; Maximum values to report different on
   maxNumDataToPrint = 6

   ;; Get types of the values
   baseType = Size(baseVal, /TNAME)
   compType = Size(compVal, /TNAME)
   
   ;; Ensure that the types of the two values are the same type
   IF (baseType NE compType) THEN BEGIN
       Report_Append, report, name + ' variables have different types, ' + $ 
         'baseline type = "' + baseType + '", ' + $
         'comparison type = "' + compType
       RETURN, 1
   ENDIF

   ;; Make sure the two types have the same type
   baseNumElems = Size(baseVal, /N_ELEMENTS)
   compNumElems = Size(compVal, /N_ELEMENTS)
   IF (baseNumElems NE compNumElems) THEN BEGIN
       Report_Append, report, name + ' variables have different sizes, ' + $ 
         'baseline size = ' + strtrim(string(baseNumElems), 1) + ', ' + $
         'comparison size = ' + strtrim(string(compNumElems), 1)
       RETURN, 1
   ENDIF

   SWITCH baseType OF
       'STRING':
       'BYTE':
       'INT':
       'LONG':
       'UINT':
       'ULONG':
       'LONG64':
       'FLOAT':
       'DOUBLE': BEGIN
           ;; Compare differently for arrays than scalars
           IF (baseNumElems GT 1) THEN BEGIN
               ;; Compare array number and string types

               SWITCH baseType OF
                   'FLOAT':
                   'DOUBLE': BEGIN
                       diffCount = Relative_Compare( baseVal, compVal, $
                                                     WHERE_DIFF=whereDiff, $
                                                     RELATIVE_DIFF=relDiff, $
                                                     _EXTRA=extra )
                       BREAK
                   END
                   ELSE: BEGIN
                       whereDiff = Where(baseVal NE compVal, diffCount)
                   END
               ENDSWITCH

               IF (diffCount GT 0) THEN BEGIN
                   ;; Output some of the values for diagnostics reasons
                   numDiffToPrint = maxNumDataToPrint < diffCount
                   IF N_Elements(relDiff) GT 0 THEN BEGIN
                       sortedWhereDiff = whereDiff[ Reverse(Sort( relDiff[whereDiff] )) ]
                       someDiffs = sortedWhereDiff[0:numDiffToPrint - 1]
                   ENDIF ELSE BEGIN
                       someDiffs = whereDiff[0:numDiffToPrint - 1]
                   ENDELSE

                   Report_Append, report, name + ' arrays of type "' + $
                     baseType + '" are different'

                   Report_Append, report, 'Out of a total of ' + $
                     strtrim(string(baseNumElems), 1) + $
                     ' values in arrays, showing ' + $
                     strtrim(string(numDiffToPrint), 1) + $
                     ' of ' + strtrim(string(diffCount), 1) + $
                     ' different values:'
               
                   Report_Append, report, 'indexes of conflicts:'
                   Report_Append, report, String(someDiffs, /PRINT)
                   
                   Report_Append, report, 'baseline values:'
                   Report_Append, report, String(baseVal[someDiffs], /PRINT)
                   
                   Report_Append, report, 'comparison values:'
                   Report_Append, report, String(compVal[someDiffs], /PRINT)

                   IF N_Elements(relDiff) GT 0 THEN BEGIN
                       Report_Append, report, 'relative differences:'
                       Report_Append, report, String(relDiff[someDiffs], /PRINT)
                   ENDIF

                   RETURN, diffCount
               ENDIF ;; end of if diffCount GT 0
               
           ENDIF ELSE BEGIN
               ;; Compare scalar number and string types
               
               SWITCH baseType OF
                   'FLOAT':
                   'DOUBLE': BEGIN
                       hasDifferences = Relative_Compare( baseVal, compVal, $
                                                          RELATIVE_DIFF=relDiff, $
                                                          _EXTRA=extra )
                       BREAK
                   END
                   ELSE: BEGIN
                       hasDifferences = baseVal NE compVal
                   END
               ENDSWITCH

               IF (hasDifferences) THEN BEGIN
                   Report_Append, report, name + ' values of type "' + $
                     baseType + '" are different. '
                   Report_Append, report, 'baseline value:'
                   Report_Append, report, String(baseVal, /PRINT)
                   Report_Append, report, 'comparison value:'
                   Report_Append, report, String(compVal, /PRINT)

                   IF N_Elements(relDiff) GT 0 THEN BEGIN
                       Report_Append, report, 'relative difference:'
                       Report_Append, report, String(relDiff, /PRINT)
                   ENDIF

                   RETURN, 1
               ENDIF

           ENDELSE

           BREAK
       END
       'COMPLEX': BEGIN
           cRealResult = Compare_Variables( Float(baseVal), $
                                            Float(compVal), $
                                            name + ":REAL", checkAll, $
                                            REPORT=cmplxRealReport, $
                                            _EXTRA=extra )
           cImagResult = Compare_Variables( Imaginary(baseVal), $
                                            Imaginary(compVal), $
                                            name + ":IMAG", checkAll, $
                                            REPORT=cmplxImagReport, $
                                            _EXTRA=extra )

           IF(cRealResult GT 0 OR cImagResult GT 0) THEN BEGIN
               Report_Append, report, cmplxRealReport
               Report_Append, report, cmplxImagReport
               
               RETURN, 1
           END

           BREAK
       END
       'DCOMPLEX': BEGIN
           dcRealResult = Compare_Variables( Double(baseVal), $
                                             Double(compVal), $
                                             name + ":REAL", checkAll, $
                                             REPORT=cmplxRealReport, $
                                             _EXTRA=extra )
           dcImagResult = Compare_Variables( Imaginary(baseVal), $
                                             Imaginary(compVal), $
                                             name + ":IMAG", checkAll, $
                                             REPORT=cmplxImagReport, $
                                             _EXTRA=extra )

           IF(dcRealResult GT 0 OR dcImagResult GT 0) THEN BEGIN
               Report_Append, report, cmplxRealReport
               Report_Append, report, cmplxImagReport
               
               RETURN, 1
           END

           BREAK
       END
       'STRUCT': BEGIN
           RETURN, Compare_Structs( baseVal, compVal, name, $
                                    checkAll, REPORT=report, $
                                    _EXTRA=extra )
           BREAK
       END
       'POINTER': BEGIN
           IF (Ptr_Valid(baseVal) AND Ptr_Valid(compVal)) THEN BEGIN
               ;; Dereference pointers and compare
               RETURN, Compare_Variables( *baseVal, *compVal, $
                                          name + ':POINTER', checkAll, $
                                          REPORT=report, $
                                          _EXTRA=extra )
           ENDIF ELSE BEGIN
               IF (Ptr_Valid(baseVal) AND NOT Ptr_Valid(compVal)) THEN BEGIN
                   Report_Append, report, name + ' has an invalid ' + $
                     'comparison pointer while the baseline pointer is valid'
                   
                   RETURN, 1
               ENDIF
               IF (NOT Ptr_Valid(baseVal) AND Ptr_Valid(compVal)) THEN BEGIN
                   Report_Append, report, name + ' has an invalid ' + $
                     'baseline pointer while the comparison pointer is valid'
                   
                   RETURN, 1
               ENDIF              
           ENDELSE

           ;; Both pointers are invalid so they are the same if we
           ;; reach here
           RETURN, 0
           
           BREAK
       END
       'OBJREF': BEGIN
           Report_Append, report, 'Recursive object compare of "' + name + $
             '" not supported'
           
           RETURN, 1
           BREAK
       END
       'UNDEFINED':
       ELSE: BEGIN
           Report_Append, report, name + ' is of unknown type "' + $
             baseType + '", ' + 'can not compare'
           
           RETURN, 1
       END
   ENDSWITCH

   ;; No differences found
   RETURN, 0
END
                          
