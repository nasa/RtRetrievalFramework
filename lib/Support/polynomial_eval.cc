#include "polynomial_eval.h"

using namespace FullPhysics;
using namespace blitz;

double Poly1d::operator()(double Value) const
{
  double res = 0;
  Array<int,1> eidx(eval_indexes());
  for(int coeff_idx = 0; coeff_idx < coeffs_.rows(); coeff_idx++) {
    res = coeffs_(eidx(coeff_idx)).value() + res * Value;
  }

  return res;
}

AutoDerivative<double> Poly1d::operator()(const AutoDerivative<double>& Value) const
{

  AutoDerivative<double> res(0);
  Array<int,1> eidx(eval_indexes());
  for(int coeff_idx = 0; coeff_idx < coeffs_.rows(); coeff_idx++) {
    res = coeffs_(eidx(coeff_idx)) + res * Value;
  }

  return res;
}

ArrayAd<double, 1> Poly1d::operator()(const ArrayAd<double, 1>& Arr) const
{
  ArrayAd<double, 1> res(Arr.rows(), Arr.number_variable());
  for(int arr_idx = 0; arr_idx < res.rows(); arr_idx++)
    res(arr_idx) = (*this)(Arr(arr_idx)); 
  return res;
}

Array<double,1> Poly1d::operator()(const blitz::Array<double,1>& Arr) const
{
  Array<double, 1> res(Arr.rows());
  for(int arr_idx = 0; arr_idx < res.rows(); arr_idx++)
    res(arr_idx) = (*this)(Arr(arr_idx)); 
  return res;
}

blitz::Array<int,1> Poly1d::eval_indexes() const
{
  Array<int, 1> res(coeffs_.rows());
  for(int idx = 0; idx < coeffs_.rows(); idx++) {
    if(decreasing_order_)
      res(idx) = idx; 
    else
      res(idx) = coeffs_.rows() - idx - 1;
  }
  return res;
}

void Poly1d::print(std::ostream& Os) const
{
  Os << "Poly1d : Polynomial: ";
  int ncoeffs = coeffs_.rows();
  Array<int,1> eidx(eval_indexes());
  for(int coeff_idx = 0; coeff_idx < ncoeffs; coeff_idx++) {
    double coeff = coeffs_(eidx(coeff_idx)).value();
    int power = ncoeffs - coeff_idx - 1;
    if(coeff_idx > 0) {
      coeff < 0.0 ? Os << " - " : Os << " + ";
      Os << abs(coeff);
    } else {
      Os << coeff;
    }
    if(power > 0) Os << "x";
    if(power > 1) Os << "^" << power;
  }
}
