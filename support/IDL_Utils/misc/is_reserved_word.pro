function is_reserved_word, word_str

   RESERVED_WORDS = [ 'AND', 'BEGIN', 'BREAK', $
                      'CASE', 'COMMON', 'COMPILE_OPT', $
                      'CONTINUE', 'DO', 'ELSE', $
                      'END', 'ENDCASE', 'ENDELSE', $
                      'ENDFOR', 'ENDIF', 'ENDREP', $
                      'ENDSWITCH', 'ENDWHILE', 'EQ', $
                      'FOR', 'FORWARD_FUNCTION', 'FUNCTION', $
                      'GE', 'GOTO', 'GT', $
                      'IF', 'INHERITS', 'LE', $
                      'LT', 'MOD', 'NE', $
                      'NOT', 'OF', 'ON_IOERROR', $
                      'OR', 'PRO', 'REPEAT', $
                      'SWITCH', 'THEN', 'UNTIL', $
                      'WHILE', 'XOR' ]
   
   where_match = where(strupcase(word_str) eq RESERVED_WORDS)
   if where_match[0] ne -1 then $
     return, 1B $
   else $
     return, 0B

end
