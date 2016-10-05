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
;;    Get_Reg_Env_Params
;;
;; FILE NAME
;;    get_reg_env_params.pro
;;
;; ABSTRACT
;;    Retrieves Regression Parameters from Enviromental Variables
;;       
;; AUTHOR
;;    James McDuffie

PRO Get_Reg_Env_Params, $
                        BASELINE_FILENAME    = baselineFilename, $
                        COMPARISON_FILENAME  = comparisonFilename, $
                        BASELINE_DIRECTORY   = baselineDirectory, $
                        COMPARISON_DIRECTORY = comparisonDirectory, $
                        TESTCASE_DIRECTORY   = testcaseDirectory, $
                        CHECK_ALL            = checkAll, $
                        INSPECT              = inspect

   baselineFilename    = $
     GetEnv('REG_BASELINE_FILENAME') ? GetEnv('REG_BASELINE_FILENAME') : ''
     
   comparisonFilename  = $
     GetEnv('REG_COMPARISON_FILENAME') ? GetEnv('REG_COMPARISON_FILENAME') : ''

   baselineDirectory   = $
     GetEnv('REG_BASELINE_DIRECTORY') ? GetEnv('REG_BASELINE_DIRECTORY') : ''

   comparisonDirectory = $
     GetEnv('REG_COMPARISON_DIRECTORY') ? GetEnv('REG_COMPARISON_DIRECTORY') : ''

   testcaseDirectory = $
     GetEnv('TESTCASE_WORKING_DIRECTORY') ? GetEnv('TESTCASE_WORKING_DIRECTORY') : ''

   checkAll            = $
     GetEnv('REG_CHECK_ALL') ? 1B : 0B

   inspect             = $
     GetEnv('REG_INSPECT') ? 1B : 0B

END
