      pro getparam, unitw, header, ncol, fname, param, nfile
;
;  Procedure to prompt the user for the names of the parameters to be read
;  and then return the selected parameter names in the array PARAM.
;
;  Inputs:
;     unitw            Logical Unit Number for writing log file.
;     nfile            Number of input files
;     fname(nfile)     Names of the input files
;     ncol(nfile)      Number of columns in each input file
;     header(nfile,*)  Array of column labels/parameters for each input file
;
;  Output:
;     param(nfile)     Names of selected parameters
;
;  Note that if the user inputs a "#" as the first character of the selected
;  parameter, the procedure will adopt this parameter name for all subsequent
;  input files, requiring no further user inputs.
;
      sss=string(" ")
      flag=0
      for ifile=0,nfile-1 do begin
         if flag eq 0  then begin
            for j=0,ncol(ifile),6 do print,format='(6a13)',header(ifile,j:j+5)
            print,format='($,"Enter parameter from  ",a)',fname(ifile)
            read,sss
            printf,unitw,sss
            sss=strtrim(sss,2)
            if strmid(sss, 0, 1) eq "#"  then  flag=1
         endif 
         param(ifile)=strmid(sss,flag,99)
      endfor
      return
      end

