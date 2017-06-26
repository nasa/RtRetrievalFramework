// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "linear_algebra.h"
%}

namespace FullPhysics {

blitz::Array<double, 1> solve_least_squares(const blitz::Array<double, 2>& A, 
                                            const blitz::Array<double, 1>& B, 
                                            double Rcond = 1e-12);
void svd(const blitz::Array<double, 2>& A, blitz::Array<double, 1>& S,
         blitz::Array<double, 2>& U, blitz::Array<double, 2>& VT);
blitz::Array<double, 2> generalized_inverse(const blitz::Array<double, 2>& A, 
                    double Rcond = std::numeric_limits<double>::epsilon());
blitz::Array<double, 1> solve_least_squares_qr(const blitz::Array<double, 2>& A, 
                                               const blitz::Array<double, 1>& B);

}
