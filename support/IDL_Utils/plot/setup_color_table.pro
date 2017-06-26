FUNCTION getBlackIndex
  Return, 1
END

FUNCTION getRedIndex
  Return, 2
END


PRO Setup_Color_Table
;  0  = [255, 255, 255]    white
;  1  = [  0,   0,   0]    black 
;  2  = [255,   0,   0]    red   
;  3  = [255, 128,   0]    orange
;  4  = [128, 128,   0]    gold
;  5  = [  0, 255,   0]    bright green
;  6  = [  0, 128, 128]    blue green
;  7  = [  0,   0, 255]    blue
;  8  = [128,   0, 255]    purple
;  9  = [255,   0, 255]    pink
;  10 = [128,   0, 128]    brown
;  11 = [  0, 255, 255]    cyan
;  12 = [255, 255,   0]    yellow

   red =   [255, 0, 255, 255, 128,   0,   0,   0, 128, 255, 128,   0, 255]
   green = [255, 0,   0, 128, 128, 255, 128,   0,   0,   0,   0, 255, 255]
   blue =  [255, 0,   0,   0,   0,   0, 128, 255, 255, 255, 128, 255,   0]

   num_colors = n_elements(red)
   dum = bytarr(256 - num_colors)

   red   = [byte(red),   dum]
   green = [byte(green), dum]
   blue  = [byte(blue),  dum]

   tvlct, red, green, blue
END
