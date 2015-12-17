#ifndef RF_GAULEG_H
#define RF_GAULEG_H

#include <blitz/array.h>
#include "auto_derivative.h"

namespace FullPhysics {

  void rf_gauleg(AutoDerivative<double> X1, AutoDerivative<double> X2, blitz::Array<AutoDerivative<double>, 1>& X, blitz::Array<AutoDerivative<double>, 1>& W, int KNUM);

}

#endif
