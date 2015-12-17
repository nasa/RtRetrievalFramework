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
;;    Relative_Compare
;;
;; FILE NAME
;;    relative_compare.pro
;;
;; ABSTRACT
;;    Performs a fractional comparison using and returns true if the
;;    values are indeed different and false otherwise.
;;       
;; AUTHOR
;;    James McDuffie

FUNCTION Relative_Compare, value1, value2, $
                           WHERE_DIFF=whereDiff, $
                           RELATIVE_DIFF=relDiff, $
                           TOLERANCE=tolerance, $
                           RMS=useRMS


   ;; Constants
   IF ( Size(value1, /TNAME) EQ 'DOUBLE' ) THEN BEGIN
       eps = (MachAr(/DOUBLE)).eps
       xmin = (MachAr(/DOUBLE)).xmin
   ENDIF ELSE BEGIN
       eps = (MachAr()).eps
       xmin = (MachAr()).xmin
   ENDELSE

   IF (N_Elements(tolerance) GT 0) THEN BEGIN
      eps = tolerance
      print, eps, " ", tolerance
  ENDIF

   ;; Add smallest value possible so we do not divide by zero
   mValue1 = value1 + 2*xmin
   mValue2 = value2 + 2*xmin

   absDiff  = abs( mValue1 - mValue2 )

   IF Keyword_Set(useRMS) THEN BEGIN
       relDiff = absDiff / Sqrt( Mean( (mValue1)^2 ) )
   ENDIF ELSE BEGIN
       relDiff = absDiff / abs(mValue1)
   ENDELSE

   whereDiff = Where(relDiff GT eps, diffCount)

   RETURN, diffCount

END
