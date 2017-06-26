      subroutine gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)
 
!           this set of routines computes the nodes t(j) and weights
!        w(j) for gaussian-type quadrature rules with pre-assigned
!        nodes.  these are used when one wishes to approximate
!
!                 integral (from a to b)  f(x) w(x) dx
!
!                              n
!        by                   sum w  f(t )
!                             j=1  j    j
!
!        (note w(x) and w(j) have no connection with each other.)
!        here w(x) is one of six possible non-negative weight
!        functions (listed below), and f(x) is the
!        function to be integrated.  gaussian quadrature is particularly
!        useful on infinite intervals (with appropriate weight
!        functions), since then other techniques often fail.
!
!           ASSOCIATED with each weight function w(x) is a set of
!        orthogonal polynomials.  the nodes t(j) are just the zeroes
!        of the proper n-th degree polynomial.
!
!     input parameters (all real numbers are in double precision)
!
!        kind     an integer between 1 and 6 giving the type of
!                 quadrature rule:
!
!        kind = 1:  legendre quadrature, w(x) = 1 on (-1, 1)
!        kind = 2:  chebyshev quadrature of the first kind
!                   w(x) = 1/sqrt(1 - x*x) on (-1, +1)
!        kind = 3:  chebyshev quadrature of the second kind
!                   w(x) = sqrt(1 - x*x) on (-1, 1)
!        kind = 4:  hermite quadrature, w(x) = EXP(-x*x) on
!                   (-infinity, +infinity)
!        kind = 5:  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**
!                   beta on (-1, 1), alpha, beta .gt. -1.
!                   note: kind=2 and 3 are a special case of this.
!        kind = 6:  generalized laguerre quadrature, w(x) = EXP(-x)*
!                   x**alpha on (0, +infinity), alpha .gt. -1
!
!        n        the number of points used for the quadrature rule
!        alpha    real parameter used only for gauss-jacobi and gauss-
!                 laguerre quadrature (otherwise use 0.d0).
!        beta     real parameter used only for gauss-jacobi quadrature--
!                 (otherwise use 0.d0)
!        kpts     (integer) normally 0, unless the left or right end-
!                 point (or both) of the interval is required to be a
!                 node (this is called gauss-radau or gauss-lobatto
!                 quadrature).  then kpts is the number of fixed
!                 endpoints (1 or 2).
!        endpts   real array of length 2.  contains the values of
!                 any fixed endpoints, if kpts = 1 or 2.
!        b        real scratch array of length n
!
!     output parameters (both double precision arrays of length n)
!
!        t        will contain the desired nodes.
!        w        will contain the desired weights w(j).
!
!        original version 20 jan 1975 from stanford
!        modified 21 dec 1983 by eric grosse
!          imtql2 => gausq2
!          hex constant => d1mach (from core library)
!          compute pi using datan
!          removed accuracy claims, description of method
!          added single precision version
!
      use linear_algebra_module
      double precision b(n), t(n), w(n), endpts(2), muzero, t1,&
       gam, solve, DSQRT, alpha, beta
!
      call class (kind, n, alpha, beta, b, t, muzero)
!
!           the matrix of coefficients is assumed to be symmetric.
!           the array t contains the diagonal elements, the array
!           b the off-diagonal elements.
!           make appropriate changes in the lower right 2 by 2
!           submatrix.
!
      if (kpts.eq.0)  go to 100
      if (kpts.eq.2)  go to  50
!
!           if kpts=1, only t(n) must be changed
!
      t(n) = solve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
      go to 100
!
!           if kpts=2, t(n) and b(n-1) must be recomputed
!
   50 gam = solve(endpts(1), n, t, b)
      t1 = ((endpts(1) - endpts(2))/(solve(endpts(2), n, t, b) - gam))
      b(n-1) = DSQRT(t1)
      t(n) = endpts(1) + gam*t1
!
!           note that the indices of the elements of b run from 1 to n-1
!           and thus the value of b(n) is arbitrary.
!           now compute the eigenvalues of the symmetric tridiagonal
!           matrix, which has been modified as necessary.
!           the method used is a ql-type method with origin shifting
!
  100 w(1) = 1.0d0
      do 105 i = 2, n
  105    w(i) = 0.0d0
!
      call gausq2 (n, t, b, w, ierr)
      do 110 i = 1, n
  110    w(i) = muzero * w(i) * w(i)
!
      return
      end subroutine gaussq
!
!
!
      double precision function solve(shift, n, a, b)
!
!       this procedure performs elimination to solve for the
!       n-th component of the solution delta to the equation
!
!             (jn - shift*identity) * delta  = en,
!
!       where en is the vector of all zeroes except for 1 in
!       the n-th position.
!
!       the matrix jn is symmetric tridiagonal, with diagonal
!       elements a(i), off-diagonal elements b(i).  this equation
!       must be solved to obtain the appropriate changes in the lower
!       2 by 2 submatrix of coefficients for orthogonal polynomials.
!
!
      double precision shift, a(n), b(n), alpha
!
      alpha = a(1) - shift
      nm1 = n - 1
      do 10 i = 2, nm1
   10    alpha = a(i) - shift - b(i-1)**2/alpha
      solve = 1.0d0/alpha
      return
      end 
!
!
!
      subroutine class(kind, n, alpha, beta, b, a, muzero)
!
!           this procedure supplies the coefficients a(j), b(j) of the
!        recurrence relation
!
!             b p (x) = (x - a ) p   (x) - b   p   (x)
!              j j            j   j-1       j-1 j-2
!
!        for the various classical (normalized) orthogonal polynomials,
!        and the zero-th moment
!
!             muzero = integral w(x) dx
!
!        of the given polynomial's weight function w(x).  since the
!        polynomials are orthonormalized, the tridiagonal matrix is
!        guaranteed to be symmetric.
!
!           the input parameter alpha is used only for laguerre and
!        jacobi polynomials, and the parameter beta is used only for
!        jacobi polynomials.  the laguerre and jacobi polynomials
!        require the gamma function.
!

      use linear_algebra_module
      double precision a(n), b(n), muzero, alpha, beta
!     double precision abi, a2b2, dgamma, pi, DSQRT, ab
      double precision abi, a2b2, pi, DSQRT, ab
!
      pi = 4.0d0 * DATAN(1.0d0)
      nm1 = n - 1
      go to (10, 20, 30, 40, 50, 60), kind
!
!              kind = 1:  legendre polynomials p(x)
!              on (-1, +1), w(x) = 1.
!
   10 muzero = 2.0d0
      do 11 i = 1, nm1
         a(i) = 0.0d0
         abi = i
   11    b(i) = abi/DSQRT(4*abi*abi - 1.0d0)
      a(n) = 0.0d0
      return
!
!              kind = 2:  chebyshev polynomials of the first kind t(x)
!              on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
!
   20 muzero = pi
      do 21 i = 1, nm1
         a(i) = 0.0d0
   21    b(i) = 0.5d0
      b(1) = DSQRT(0.5d0)
      a(n) = 0.0d0
      return
!
!              kind = 3:  chebyshev polynomials of the second kind u(x)
!              on (-1, +1), w(x) = sqrt(1 - x*x)
!
   30 muzero = pi/2.0d0
      do 31 i = 1, nm1
         a(i) = 0.0d0
   31    b(i) = 0.5d0
      a(n) = 0.0d0
      return
!
!              kind = 4:  hermite polynomials h(x) on (-infinity,
!              +infinity), w(x) = EXP(-x**2)
!
   40 muzero = DSQRT(pi)
      do 41 i = 1, nm1
         a(i) = 0.0d0
   41    b(i) = DSQRT(i/2.0d0)
      a(n) = 0.0d0
      return
!
!              kind = 5:  jacobi polynomials p(alpha, beta)(x) on
!              (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
!              beta greater than -1
!
   50 ab = alpha + beta
      abi = 2.0d0 + ab
      muzero = 2.0d0 ** (ab + 1.0d0) * dgamma(alpha + 1.0d0) * dgamma(&
       beta + 1.0d0) / dgamma(abi)
      a(1) = (beta - alpha)/abi
      b(1) = DSQRT(4.0d0*(1.0d0 + alpha)*(1.0d0 + beta)/((abi + 1.0d0)*&
        abi*abi))
      a2b2 = beta*beta - alpha*alpha
      do 51 i = 2, nm1
         abi = 2.0d0*i + ab
         a(i) = a2b2/((abi - 2.0d0)*abi)
   51    b(i) = DSQRT (4.0d0*i*(i + alpha)*(i + beta)*(i + ab)/&
         ((abi*abi - 1)*abi*abi))
      abi = 2.0d0*n + ab
      a(n) = a2b2/((abi - 2.0d0)*abi)
      return
!
!              kind = 6:  laguerre polynomials l(alpha)(x) on
!              (0, +infinity), w(x) = EXP(-x) * x**alpha, alpha greater
!              than -1.
!
   60 muzero = dgamma(alpha + 1.0d0)
      do 61 i = 1, nm1
         a(i) = 2.0d0*i - 1.0d0 + alpha
   61    b(i) = DSQRT(i*(i + alpha))
      a(n) = 2.0d0*n - 1 + alpha
      return
      end subroutine class
!
!
      subroutine gausq2(n, d, e, z, ierr)
!
!     this subroutine is a translation of an algol procedure,
!     num. math. 12, 377-383(1968) by martin and wilkinson,
!     as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
!     this is a modified version of the 'eispack' routine imtql2.
!
!     this subroutine finds the eigenvalues and first components of the
!     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
!     method.
!
!     on input:
!
!        n is the order of the matrix;
!
!        d contains the diagonal elements of the input matrix;
!
!        e contains the subdiagonal elements of the input matrix
!          in its first n-1 positions.  e(n) is arbitrary;
!
!        z contains the first row of the identity matrix.
!
!      on output:
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1, 2, ..., ierr-1;
!
!        e has been destroyed;
!
!        z contains the first components of the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.  if an error exit is
!          made, z contains the eigenvectors ASSOCIATED with the stored
!          eigenvalues;
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     ------------------------------------------------------------------
!

      use linear_algebra_module
      integer i, j, k, l, m, n, ii, mml, ierr
      real*8 d(n), e(n), z(n), b, c, f, g, p, r, s, machep
!     real*8 DSQRT, DABS, dsign, d1mach
      real*8 DSQRT, DABS, DSIGN
!
      machep=D1MACH(4)
!
      ierr = 0
      if (n .eq. 1) go to 1001
!
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
!     :::::::::: look for small sub-diagonal element ::::::::::
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            if (DABS(e(m)) .le. machep * (DABS(d(m)) + DABS(d(m+1))))&
               go to 120
  110    continue
!
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
!     :::::::::: form shift ::::::::::
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = DSQRT(g*g+1.0d0)
         g = d(m) - p + e(l) / (g + DSIGN(r, g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
!
!     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (DABS(f) .lt. DABS(g)) go to 150
            c = g / f
            r = DSQRT(c*c+1.0d0)
            e(i+1) = f * r
            s = 1.0d0 / r
            c = c * s
            go to 160
  150       s = f / g
            r = DSQRT(s*s+1.0d0)
            e(i+1) = g * r
            c = 1.0d0 / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
!     :::::::::: form first component of vector ::::::::::
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
  200       z(i) = c * z(i) - s * f
!
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
  240 continue
!
!     :::::::::: order eigenvalues and eigenvectors ::::::::::
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
!
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
!
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         p = z(i)
         z(i) = z(k)
         z(k) = p
  300 continue
!
      go to 1001
!     :::::::::: set error -- no convergence to an
!                eigenvalue after 30 iterations ::::::::::
 1000 ierr = l
 1001 return
!     :::::::::: last card of gausq2 ::::::::::
      end subroutine gausq2
