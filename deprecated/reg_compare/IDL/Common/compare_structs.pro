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
;;    Compare_Structs
;;
;; FILE NAME
;;    compare_structs.pro
;;
;; ABSTRACT
;;    Compares two IDL structures
;; 
;; AUTHOR
;;    James McDuffie
;;
;; DESCRIPTION
;;    Compares two IDL structures without knowing anything about how
;;    they are formatted. It returns true if the structs have differences
;;    and false otherwise. Additionally, a textual report can be
;;    retrieved through a named parameter. 
;;    The optional TAGS_TO_CHECK keyword is used whenever the user
;;    does not want all tags check

FUNCTION Compare_Structs, $
                          s_baseStruct,  $
                          s_compStruct,  $
                          name,          $
                          checkAll,      $
                          REPORT=report, $
                          TAGS_TO_CHECK=tagsToCheck, $
                          _EXTRA=extra


   ;; Set the report to state that the structures are equal this will be
   ;; overwritten if the function detects that they are indeed not
   report = ''
   numDifferences = 0

   ;; Initialize the list of tags that will be used
   IF (N_Elements(tagsToCheck) GT 0) THEN BEGIN
       baseTags = tagsToCheck
       compTags = tagsToCheck
   ENDIF ELSE BEGIN
       baseTags = tag_names(s_baseStruct)
       compTags = tag_names(s_compStruct)
   ENDELSE

   ;; Compare number of variables within the structures
   IF(N_Elements(baseTags) NE N_Elements(compTags)) THEN BEGIN
       ;; Log what has happened and return failure
       report = 'Number of variables in struct "' + name + '" are not equal'
       RETURN, 1
   ENDIF

   ;; Check each pair of tags
   FOR i = 0, N_Elements(baseTags) - 1 DO BEGIN

       ;; Make sure names of variables are equal
       IF(baseTags[i] NE compTags[i]) THEN BEGIN
           report = 'Variable names of "' + name + '" structs at index ' + $
             strtrim(string(i), 1) + ' are not the same.' + $ 
             ' baseline name: ' +  baseTags[i] + $
             ' comparison name: ' + compTags[i]

           RETURN, 1         
       ENDIF

       tagIndex = (Where(tag_names(s_baseStruct) EQ baseTags[i]))[0]

       IF (tagIndex EQ -1) THEN BEGIN
           report = 'Error finding index of tag ' + baseTags[i]
           RETURN, 1
       ENDIF

       numVarDiffs = $
         Compare_Variables( s_baseStruct.(tagIndex), s_CompStruct.(tagIndex), $
                            name + '.' + baseTags[i], $
                            checkAll, REPORT=varReport, $
                            _EXTRA=extra )
       
       IF(numVarDiffs GT 0) THEN BEGIN
           Report_Append, report, varReport
           numDifferences = numDifferences + 1
           IF(NOT checkAll) THEN RETURN, numDifferences
       ENDIF
   ENDFOR

   ;; Return number of differences in a check all situation or 0 if no
   ;; differences
   RETURN, numDifferences

END

