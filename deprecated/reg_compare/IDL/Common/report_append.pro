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
;;    Report_Append
;;
;; FILE NAME
;;    report_append.pro
;;
;; ABSTRACT
;;    Convenience method to concatenate stringsToAdd to report and add
;;    a newline.
;; 
;; AUTHOR
;;    James McDuffie

PRO Report_Append, report, stringsToAdd

   ;; Get a string representing a newline
   newLine = string(byte(10))

   FOR i = 0, N_Elements(stringsToAdd) - 1 DO BEGIN
       report = report + stringsToAdd[i] + newLine
   ENDFOR

END
