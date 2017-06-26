;+
;;##################### TES Level 1B IDL Code ########################
;;
;; ABSTRACT
;;    Compares two OCO L2 Matrix files
;;       
;; AUTHOR
;;    James McDuffie

FUNCTION Compare_Matrix_Files, $
  baselineFilename, $
  comparisonFilename, $
  checkAll, $
  REPORT=report

   o_baseLineObj = Obj_New('MATRIX_FILE', baselineFilename)
   o_comparisonObj = Obj_New('MATRIX_FILE', comparisonFilename)

   ;; No differences so far
   numDifferences = 0
   report = ''

   objectMethods = [ $
                   'Get_Num_Rows', $
                   'Get_Num_Columns', $
                   'Get_Header_Keyword_Names', $
                   'Get_File_ID', $
                   'Get_File_Type', $
                   'Get_Data' $
                   ]

   ;; Specify methods to check the size of array data types
   sizeCheckMethods = StrArr( N_Elements(objectMethods) )
   sizeCheckMethods[Where(objectMethods EQ 'Get_Data')] = 'Get_Num_Rows'

   ;; Use RMS for spectra type comparisons
   rmsMethods = BytArr( N_Elements(objectMethods) )
   rmsMethods[Where(objectMethods EQ 'Get_Data')] = 1B

   numDiffs =  Compare_Objects( o_baselineObj, $
                                o_comparisonObj, $
                                objectMethods, $
                                checkAll, $
                                SIZE_CHECK_METHODS=sizeCheckMethods, $
                                RMS_COMPARE_METHODS=rmsMethods, $
                                REPORT=report )

   Return, numDiffs
   
END
