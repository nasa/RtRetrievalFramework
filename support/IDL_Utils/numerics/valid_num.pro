;+
; NAME: 
;     VALID_NUM
; PURPOSE:               
;     Check if a string is a valid number representation.
; EXPLANATION:              
;     The input string is parsed for characters that may possibly
;     form a valid number.  It is more robust than simply checking
;     for an IDL conversion error because that allows strings such
;     as '22.3qwert' to be returned as the valid number 22.3
;     See also the original NUM_CHK which returns the status in 
;     the opposite sense.
;
; CALLING SEQUENCE: 
;     IDL> status = valid_num(string  [,value]  [,/integer])
;    
; Inputs      : string  -  the string to be tested
;               
; Opt. Inputs : None
;               
; Outputs     : The function returns 1 for valid, 0 for invalid number
;               
; Opt. Outputs: value	- The value the string decodes to.  This will be
;			  returned as a double precision number unless /INTEGER
;			  is present, in which case a long integer is returned.
;               
; Keywords    : Integer   -  if present code checks specifically for an integer.
;
; Calls       : None
;               
; Restrictions: None
;               
; Category    : Utilities, Numerical
;               
; Prev. Hist. : Small changes from NUM_CHK by Andrew Bowen, 
;                                             Tessella Support Services, 8/3/93
;
; Written     : CDS version by C D Pike, RAL, 24-May-93
;               
; Modified    : Version 1, C D Pike, RAL, 24-May-93
;		Version 2, William Thompson, GSFC, 14 October 1994
;			Added optional output parameter VALUE to allow
;			VALID_NUM to replace STRNUMBER in FITS routines.
;
; Version     : Version 1  24-May-93
;	Converted to IDL V5.0   W. Landsman   September 1997
;-            

FUNCTION valid_num, string, value, INTEGER=integer

  ;; Ensure we only use a scalar
  string = string[0]

		;**** Set defaults for keyword ****
  IF NOT (KEYWORD_SET(integer)) THEN integer=0

		;**** arrays of legal characters ****
  numbers 	= '0123456789'
  signs 	= '+-'
  decimal 	= '.'
  exponents 	= 'ED'

		;**** trim leading and trailing blanks/compress white ****
		;**** space and convert any exponents to uppercase.   ****
  numstr = strupcase(strtrim(strcompress(string),2))

		;**** length of input string ****
  len = strlen(numstr)

  ok = 1

  if integer eq 0 then stage = 1 else stage = 6

  for i = 0, len-1 do begin

    char = strmid(numstr,i,1)

		;**** the parsing steps 1 to 8 are for floating   ****
		;**** point, steps 6 to 8, which test for a legal ****
		;**** exponent, can be used to check for integers ****

;**** The parsing structure is as follows.  Each character in the ****
;**** string is checked against the valid list at the current     ****
;**** stage.  If no match is found an error is reported.  When a  ****
;**** match is found the stage number is updated as indicated     ****
;**** ready for the next character.  The valid end points are     ****
;**** indicated in the diagram.					  ****
;
;Stage	1		2		3		4
;
;Valid	sign	--> 2	dec-pt	--> 3	digit	--> 5	dec-pt	--> 5
;  "	dec-pt	--> 3	digit	--> 4			digit	--> 4
;  "	digit	--> 4					exp't	--> 6
;  "							END
;
;Stage	5		6		7		8
;
;Valid	digit	--> 5	sign	--> 7	digit	--> 8	digit	-->8
;  "	exp't	--> 6	digit	--> 8			END
;  "	END
;

    CASE stage OF

      1 : begin
        if 		strpos(signs,char) ge 0 	then stage = 2 $
	else if 	decimal eq char 		then stage = 3 $
	else if 	strpos(numbers,char) ge 0 	then stage = 4 $
	else 		ok = 0
      end

      2 : begin
	if	 	decimal eq char 		then stage = 3 $
	else if 	strpos(numbers,char) ge 0 	then stage = 4 $
	else 		ok = 0
      end

      3 : begin
	if	 	strpos(numbers,char) ge 0 	then stage = 5 $
	else 		ok = 0
      end

      4 : begin
	if	 	decimal eq char 		then stage = 5 $
	else if 	strpos(numbers,char) ge 0 	then stage = 4 $
	else if		strpos(exponents,char) ge 0	then stage = 6 $
	else 		ok = 0
      end

      5 : begin
	if	 	strpos(numbers,char) ge 0 	then stage = 5 $
	else if		strpos(exponents,char) ge 0	then stage = 6 $
	else 		ok = 0
      end

      6 : begin
        if 		strpos(signs,char) ge 0 	then stage = 7 $
	else if 	strpos(numbers,char) ge 0 	then stage = 8 $
	else 		ok = 0
      end

      7 : begin
	if	 	strpos(numbers,char) ge 0 	then stage = 8 $
	else 		ok = 0
      end

      8 : begin
	if	 	strpos(numbers,char) ge 0 	then stage = 8 $
	else 		ok = 0
      end

    ENDCASE

  end

		;**** check that the string terminated legally ****
		;**** i.e in stages 4, 5 or 8                  ****
  if (stage ne 4) and (stage ne 5) and (stage ne 8) then ok = 0

		;**** If requested, then form the value. ****

  if (n_params() eq 2) and ok then begin
	if keyword_set(integer) then value = long(string) else	$
		value = double(string)
  endif

		;**** return error status to the caller ****
  RETURN, ok


END
