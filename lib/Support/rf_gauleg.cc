#include "rf_gauleg.h"

// Gauss-Legendre Integration
void FullPhysics::rf_gauleg(AutoDerivative<double> X1, AutoDerivative<double> X2, blitz::Array<AutoDerivative<double>, 1>& X, blitz::Array<AutoDerivative<double>, 1>& W, int KNUM) {

  X.resize(KNUM);
  W.resize(KNUM);

  int I2, M;
  double P1,P2,P3,PP,Z,Z1,ZM1,ZP1,HELP,ZH;
  const double EPS = 3.e-14;

  AutoDerivative<double> XL = 0.5e0 * ( X2 - X1 );
  M=(KNUM+1)/2;
  for(int I1 = 1; I1 <= M; I1++) {
    I2 = KNUM + 1 - I1;
    Z = cos(3.141592654E0 * (I1-.25E0) / (KNUM +.5E0));
    do {
      P1=1.E0;
      P2=0.E0;
      for(int L=1; L <= KNUM; L++) {
	P3=P2;
	P2=P1;
	P1=((2.E0*L-1.E0)*Z*P2-(L-1.E0)*P3)/L;
      }
      PP=KNUM*(Z*P1-P2)/(Z*Z-1.E0);
      Z1=Z;
      Z=Z1-P1/PP;
    } while  (fabs(Z-Z1) > EPS);
    ZH = 0.5e0 * Z;
    ZP1 = 0.5e0 + ZH;
    ZM1 = 0.5e0 - ZH;
    X(I1-1) = ZM1 * X2 + ZP1 * X1;
    X(I2-1) = ZP1 * X2 + ZM1 * X1;
    HELP = 2.E0/((1.E0-Z*Z)*PP*PP);
    W(I1-1) = HELP * XL;
    W(I2-1) = W(I1-1);
  }
}
