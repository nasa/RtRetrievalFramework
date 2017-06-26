;+
;;##################### TES Level 1B IDL Code ########################
;;
;; ABSTRACT
;;    Loads in a filename passed to it from the Perl module then loops
;;    over the file pairs and calls another function to determine the
;;    correct module to use for comparison.
;;       
;; AUTHOR
;;    James McDuffie

PRO OCO_L2_Matrix_Compare

   !Except = 0

   ;; Retrieve enviromental variables used by this module
   Get_Reg_Env_Params, CHECK_ALL=checkAll, INSPECT=inspect
   pairListFilename = GetEnv('REG_PAIR_LIST_FILENAME')

   ;; cd to location where script is being run from
   cd, GetEnv('TESTCASE_WORKING_DIRECTORY')

   ;; Error handling. Make sure that we can inform the caller of an error
   errStat = 0
;   IF NOT Keyword_Set(inspect) THEN Catch, errStat
   IF errStat NE 0 THEN BEGIN
       ON_ERROR, 2
       Catch, /CANCEL

       OpenW, lun, pairListFilename, /GET_LUN, ERROR=openErr
       PrintF, lun, 'ERROR ERROR'
       
       IF  (openErr NE 0) THEN $
         Message, /CONTINUE, !Error_State.msg $
       ELSE $
         Free_Lun, lun
       
       Message, !Error_state.msg
   ENDIF

   ;; Initial size to be read
   maxRec = 5000000L
   baseFilenames = StrArr(maxRec)
   compFilenames = StrArr(maxRec)

   OpenR, lun, pairListFilename, /GET_LUN, ERROR=openErr
   IF (openErr NE 0) THEN Message, !Error_State.msg, /NONAME
   
   ;; Read the records until end of line
   nRecords = 0L
   nLine = 1L 
   WHILE (eof(lun) NE 1) DO BEGIN
       ;; Read current record, jump to bad_record: on error
       on_ioerror, bad_record
       error = 1
       str = ''
       ReadF, lun, str

       parts = StrSplit(str, /extract)
       IF (N_Elements(parts) NE 2) THEN Message, 'Not 2 parts'

       baseFilenames[nRecords] = parts[0]
       compFilenames[nRecords] = parts[1]
       nRecords = nRecords + 1

       ;; Check for bad input
       bad_record:
       IF (error EQ 1 AND Keyword_Set(showBad)) THEN BEGIN
           print, 'Bad record encountered at line: ', nLine
       ENDIF
       
       nLine = nLine + 1
   ENDWHILE

   ;; Close pair filex
   Free_Lun, lun

   baseFilenames = baseFilenames[0:nRecords - 1]
   compFilenames = compFilenames[0:nRecords - 1]

   lineStr = ''
   FOR i=1, 80 DO $
     lineStr = lineStr + '-'

   ;; Write out difference counts to same file read list of filename
   ;; from. Add the count of differerences to the file map
   OpenW, lun, pairListFilename, /GET_LUN, ERROR=openErr

   numDiffFiles = 0
   FOR fileIndex = 0L, nRecords - 1 DO BEGIN;    
       numPairDiffs = Compare_Matrix_Files( baseFilenames[fileIndex], $
                                            compFilenames[fileIndex], $
                                            checkAll, $
                                            REPORT=report )

       IF numPairDiffs GT 0 THEN BEGIN
           print, lineStr
           print, 'Comparison of ' + baseFilenames[fileIndex] + $
             ' with ' + compFilenames[fileIndex]
           print, ''

           Print, report
           Print, ''
           Print, numPairDiffs, FORMAT='("Files have ", I0, " differences")'
           numDiffFiles = numDiffFiles + 1
           IF (NOT checkAll) THEN Return_Result, 1
       ENDIF

       PrintF, lun, String(baseFilenames[fileIndex], compFilenames[fileIndex], numPairDiffs, FORMAT='(A0, " ", A0, " ", I0)')
   ENDFOR

   IF (openErr NE 0) THEN $
     Message, !Error_State.msg, /NONAME $
   ELSE $
     Free_Lun, lun

   IF (numDiffFiles GT 0) THEN BEGIN
       Return_Result, 1
   ENDIF ELSE BEGIN
       Return_Result, 0
   ENDELSE
   
END
