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
;;    Compare_Objects
;;
;; FILE NAME
;;    compare_objects.pro
;;
;; ABSTRACT
;;    Compares two Objects
;; 
;; AUTHOR
;;    James McDuffie
;;
;; DESCRIPTION
;;    Compares two IDL objects by iterating through a list of method
;;    names to test. The methods named should be functions that take
;;    no arguments and return their namesake. This function returns
;;    true if the objects have differences and false otherwise.
;;    The optional size check methods keyword specifies a string array
;;    of the same size as the method names array. If an item in the
;;    size check methods array is not the empty string then it
;;    specifies the name of a method to use to check an array's
;;    size before calling the method that would return the array. This
;;    is useful when a method will try to dereference a pointer even
;;    if it is invalid.
;;    Additionally, a textual report can be retrieved through a named
;;    parameter.


FUNCTION Compare_Objects, o_baselineObj, $
                          o_comparisonObj, $
                          methodNames, $
                          checkAll, $
                          SIZE_CHECK_METHODS=sizeCheckMethods, $
                          RMS_COMPARE_METHODS=rmsCompareMethods, $
                          METHOD_NAME_PREPEND=methodNamePrepend, $
                          REPORT=report, $
                          _EXTRA=extra
  
   ;; Initially there are zero differences
   numDifferences = 0
   
   ;; Initialize report
   report = ''

   ;; Make sure that if the size check methods array is specified it
   ;; is the same size as the methods name array
   IF (N_Elements(sizeCheckMethods) GT 0) THEN BEGIN
       IF(N_Elements(sizeCheckMethods) NE N_Elements(methodNames)) THEN BEGIN
           Message, 'Size check methods keyword specified but is not the ' + $
             'same size as the method names array', /CONTINUE
           EXIT, STATUS=1
       ENDIF
   ENDIF

   ;; Make sure both objects are both valid
   IF (NOT (Obj_Valid(o_baselineObj) AND Obj_Valid(o_comparisonObj))) THEN BEGIN
       report = 'Both objects are not valid. baseline object valid: ' + $
         StrTrim(String(Obj_Valid(o_baselineObj)), 1) + $
         ', comparison object valid: ' + $
         StrTrim(String(Obj_Valid(o_comparisonObj)), 1)
       RETURN, 1
   ENDIF

   ;; Ensure that both objects are of the same class
   IF (Obj_Class(o_baselineObj) NE Obj_Class(o_comparisonObj)) THEN BEGIN
       report = 'Class type of baseline object "' + $
         Obj_Class(o_baselineObj) + '" and comparison object "' + $
         Obj_Class(o_comparisonObj) + '" are not the same'
       
       RETURN, 1
   ENDIF ELSE BEGIN
       className = Obj_Class(o_baselineObj)
   ENDELSE

   ;; Make sure that method names have been passed to the class
   numMethods = size(methodNames, /N_ELEMENTS)

   IF (numMethods EQ 0) THEN BEGIN
       report = 'No methods names passed for objects of type ' + $
         Obj_Class(o_baselineObj)
       RETURN , 1
   ENDIF

   ;; Iterate over the list of method names
   FOR i = 0, numMethods - 1 DO BEGIN
       currMethod = methodNames[i]

       baseSize = 1
       compSize = 1
       IF ( N_Elements(sizeCheckMethods) GT 0 ) THEN BEGIN
           IF ( SizeCheckMethods[i] NE '' ) THEN BEGIN

               currSizeMethod = sizeCheckMethods[i]
               baseSize = Call_Method(currSizeMethod, o_baselineObj)
               compSize = Call_Method(currSizeMethod, o_comparisonObj)
           
               IF (baseSize NE compSize) THEN BEGIN
                   Report_Append, report, 'The size methods for the data ' + $
                     'method "' + currMethod + '" are not equal. ' + $
                     'baseline data size = ' + $
                     StrTrim(String(baseSize), 1) + ', ' + $
                     'comparison data size = ' + $
                     StrTrim(String(compSize), 1)
                   IF (NOT checkAll) THEN RETURN, numDifferences
                   Continue
               ENDIF
           ENDIF
       ENDIF

       numScalDiffs = 0
       IF (baseSize GE 1 AND compSize GE 1) THEN BEGIN
           baseVal = Call_Method(currMethod, o_baselineObj)
           compVal = Call_Method(currMethod, o_comparisonObj)

           IF Keyword_Set(methodNamePrepend) THEN $
             currMethodStr = methodNamePrepend + currMethod $
           ELSE $
             currMethodStr = currMethod

           ;; See if we shall use RMS to do the comparison for the type
           useRMS = 0B
           IF( N_Elements(rmsCompareMethods) GT 0 ) THEN BEGIN
               IF( rmsCompareMethods[i] EQ 1B ) THEN BEGIN
                   useRMS = 1B
               ENDIF
           ENDIF
       
           IF (useRMS) THEN BEGIN
               numScalDiffs = Compare_Variables( baseVal, compVal, currMethodStr, $
                                                 checkAll, REPORT=sReport, $
                                                 /RMS, _EXTRA=extra )

           ENDIF ELSE BEGIN
               numScalDiffs = Compare_Variables( baseVal, compVal, currMethodStr, $
                                                 checkAll, REPORT=sReport, $
                                                 _EXTRA=extra )
           ENDELSE

       ENDIF
       
       IF (numScalDiffs GT 0) THEN BEGIN
           Report_Append, report, sReport
           numDifferences = numDifferences + numScalDiffs
           
           IF (NOT checkAll) THEN RETURN, numDifferences
       ENDIF
   ENDFOR

   RETURN, numDifferences
END

