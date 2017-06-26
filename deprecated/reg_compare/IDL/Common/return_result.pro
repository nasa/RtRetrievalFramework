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
;;    Return_Result
;;
;; FILE NAME
;;    return_result.pro
;;
;; ABSTRACT
;;    Returns the result of the regression test back to the caller
;;    through enviromental variables
;;       
;; AUTHOR
;;    James McDuffie

PRO Return_Result, resultCode

   ;; Determine if we are in debugging mode if so then do not exit
   Get_Reg_Env_Params, $
     INSPECT=inspect
  
   IF (inspect AND resultCode > 0) THEN BEGIN
       print, "Result Code = ", resultCode
       Stop
   ENDIF ELSE BEGIN
       Exit, STATUS=resultCode
   ENDELSE

END
