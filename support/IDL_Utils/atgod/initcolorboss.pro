; $Id: initcolorboss.pro,v 1.11 2008/03/14 23:02:05 fullerr Exp $
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; ColorBoss
; <p>
; This is a new set of idl routines, based on the previous
; color_manager set of routines.  However, this module has been
; brought up to date significantly.
; <p>
; Perhaps the most notable thing is the use of 24 bit color, and the
; fact that I'm dealing with loadct correctly.
; <p>
; Also note, I've finally bitten the bullet and moved to the American
; spelling of colour!
; <p>
; This new version also has the space for a transfer function for
; funny printers (e.g. pastel2)
; 
; @param defaultFile {type=String} {Optional} {default='~/idl/my_tables.tbl'}
;                    The table to get the normal color values from.
; @keyword noReset   {type=Boolean}
;                    When set, if the Color Boss has previously been
;                    set up, return without doing any work.
;
; @author Nathaniel Livesey
; @version $Revision: 1.11 $ $Date: 2008/03/14 23:02:05 $
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro InitColorBoss, defaultFile, noReset=noReset
common ColBoss,config

;; If noReset given and color boss has previously been set up, return gracefully
IF N_Elements(config) NE 0 && Keyword_Set(noReset) THEN RETURN

; If supplied, defaultFile is the color table file the user wishes to
; use by default.

if n_elements(defaultFile) eq 0 then begin
  ;; Get the location of this directory
  Help, /SOURCE, /PROCEDURES, name='initcolorboss', output=op
  
  ;; The info we want is on the last line after the spaces
  root = File_DirName((StrSplit(op[N_Elements(op)-1], ' ', /EXTRACT, count=cnt))[cnt-1], $
                      /MARK_DIRECTORY)

  defaultFile= root+"my_tables.tbl"
endif

; First ascertain whether we're in pseudo color or true color mode

if !D.name eq 'X' then begin
  Device,get_visual_name=visual
  pseudo=visual eq 'PseudoColor'
endif else pseudo=1

; Force decomposed color if not in pseudo color mode

if not pseudo then Device,/decomposed

; If there was a transfer function before, then keep hold of it.

if n_elements(config) ne 0 $
  then transferFunction=config.transferFunction $
else transferFunction=''

config={ $
  pseudo:pseudo, $              ; Flag to indicate in pseudo color mode.
  defaultCt:34, $                ; Index number for default color table.
  defaultFile:defaultFile, $    ; Filename for default ct file.
  transferFunction:transferFunction, $
  firstFree:0L}                 ; First free color (only needed in pseudo mode)

FlushCMYKMappings

end

; $Log: initcolorboss.pro,v $
; Revision 1.11  2008/03/14 23:02:05  fullerr
; Added noReset keyword
;
; Revision 1.10  2005/07/14 01:05:54  livesey
; Changed default color table
;
; Revision 1.9  2003/09/19 00:28:02  fullerr
; Added documentation
;
; Revision 1.8  2003/05/05 23:50:37  livesey
; Added call to flush cmyk mappings
;
; Revision 1.7  2001/07/03 00:08:48  livesey
; Reformatted and indented etc.
;
; Revision 1.6  2000/11/06 23:48:09  livesey
; Moved home directory

; Revision 1.5  2000/08/28 20:36:14  livesey
; Replaced /gwynedd/livesey with /gwynedd/livesey unfortunately

; Revision 1.4  2000/01/05 04:55:17  livesey
; Regular commit

; Revision 1.3  1999/10/21 17:18:40  livesey
; Added transfer function stuff to colorboss

; Revision 1.2  1999/09/23 22:28:42  livesey
; Routine commit

; Revision 1.1  1999/01/18 18:31:29  livesey
; Initial revision





