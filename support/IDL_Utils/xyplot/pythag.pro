pro pythag, a, b, c, nrow 
;
;  Inputs:
;      A(NROW)   First array of values
;      B(NROW)   Second array of values
;         NROW   Number of values
;
; Outputs:
;      C(NROW)   SQRT(A^2+B^2)
;
; Calculates SQRT(a^2+b^2) without actually performing a SQRT.
; Also, no risk of overflow or underflow.

for i=long(0),nrow-1 do begin
   aa=abs(a(i))
   bb=abs(b(i))
   p = max([aa,bb])
   if (p gt 0.0e0) then begin
      r = (min([aa,bb])/p) ^ 2
      t = 4.0e0 + r
      while (t ne 4.0e0) do begin
         s = r/t
         u = 1.0e0 + 2.0e0 * s
         p = u * p
         r = (s/u)^2 * r
         t = 4.0e0 + r
      endwhile
   endif
   c(i)=p
endfor
return
end

