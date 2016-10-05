PRO Copy_Struct_Contents, srcStruct, dstStruct


   srcTags = tag_names(srcStruct)
   dstTags = tag_names(dstStruct)

   ;; Check each source tag to copy into dest
   FOR srcTagIndex = 0, N_Elements(srcTags) - 1 DO BEGIN
       dstTagIndex = (Where(dstTags EQ srcTags[srcTagIndex]))[0]

       IF (dstTagIndex NE -1) THEN BEGIN
           dstStruct.(dstTagIndex) = srcStruct.(srcTagIndex)         
       ENDIF ELSE BEGIN
           IF size(dstStruct.(dstTagIndex), /tname) eq 'STRING' THEN BEGIN
               dstStruct.(dstTagIndex) = ''
           ENDIF ELSE BEGIN
               dstStruct.(dstTagIndex) = 0
           ENDELSE
       ENDELSE

   ENDFOR

END

