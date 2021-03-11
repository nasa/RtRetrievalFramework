#ifndef ILS_GAUSSIAN_H
#define ILS_GAUSSIAN_H
#include "ils_function.h"

namespace FullPhysics {
/****************************************************************//**
  This is an ILS function that is a Gaussian.
*******************************************************************/
class IlsGaussian : public IlsFunction {
public:
//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
  IlsGaussian(double a, const std::string& Band_name, 
	      const std::string& Hdf_band_name) 
    : band_name_(Band_name), hdf_band_name_(Hdf_band_name)
  { gamma = a / sqrt(2); }
  virtual ~IlsGaussian() {}

  virtual void ils
  (const AutoDerivative<double>& wn_center,
   const blitz::Array<double, 1>& wn, ArrayAd<double, 1>& res,
   bool jac_optimization = false) const;
  virtual void print(std::ostream& os) const {os << "IlsGaussian";}
  virtual std::string band_name() const { return band_name_; }
  virtual std::string hdf_band_name() const {return hdf_band_name_;}
private:
  double gamma;
  std::string band_name_, hdf_band_name_;
};
}

//-----------------------------------------------------------------------
// For reference, here are the other functions that use to be in the
// old Fortran code. We can easily implement and test these, but there
// doesn't seem to be much point until we actually need them:
//
/*
void ils_coefs(double wn_center, Array<double,1> wn,
  const Array<double,1>& para, int ils_function,
  Array<double,1>& ils_out)
{
  int wn_length = wn.rows();
  int num_para = para.rows();
  Range allwn(0,wn_length-1);

     // constants

  const double sqrt_log_2 = 0.83255461115769769;
  const double sqrt_PI = 1.7724538509055159;

      // local variables

  Array<double,1> thres(3);
  double a;
  double c;
  double d;
  double Peak;
  double s;
  double Cbkg;
  double Cs;
  double t;

  Array<double,1> x(wn_length);
  double b;
  Array<double,1> L1(wn_length);
  Array<double,1> L2(wn_length);
  double gamma;
  int i;

    // init

  wn.reindexSelf(TinyVector<int,1>(0)); // reindex local version only

    // generate the ILS per the selected model and parameters

  switch (ils_function) {
    case KEY_ILS_LORENTZ:
       x(allwn) = wn(allwn) - wn_center;
       if (num_para > 1)
         throw Exception("Error ! Number of ILS Parameters > 1");
       a = para(0);
       ils_out(allwn) = a / (x(allwn)*x(allwn) + a*a);
       break;

    case KEY_ILS_SUM_LORENTZ:
         // if only one parameter is given, then Lorentz6 will be used

       x(allwn) = wn(allwn) - wn_center;
       if (num_para == 1) {
          a = para(0);
          b = 1.0;
       }
       else if (num_para == 2) {
          a = para(0);
          b = para(1);
       }
       else throw Exception("Error ! Number of ILS Parameters > 2");

       L1(allwn) = (a/PI) / (x(allwn)*x(allwn) + a*a);
       L2(allwn) = pow((6.0 * a / (x(allwn)*x(allwn) + 6.0*a*a)), 6.0) /
         pow((PI / a), 6.0);
       ils_out(allwn) = b * L2(allwn) + (1.0 - b) * L1(allwn);
       break;

    case KEY_ILS_GAUSS:
       x(allwn) = wn(allwn) - wn_center;
       if (num_para > 1)
          throw Exception("Error ! Number of ILS Parameters > 1");
       a = para(0);

       gamma = a/(sqrt_log_2); // convert HWHM into e-length

       ils_out(allwn) = exp(-(x(allwn)/gamma)*(x(allwn)/gamma)) / (gamma*sqrt_PI);
       break;

    case KEY_ILS_SIGMOID:
       x(allwn) = wn(allwn) - wn_center;
       Cbkg = 0.0;
       Cs = 1.0;
       if (num_para > 3)
         throw Exception("Error ! Number of ILS Parameters > 3");
       a = para(0);
       b = para(1);
       s = para(2);

       L1(allwn) = -a/2.0+b;
       L2(allwn) = a/2.0+b;
       ils_out(allwn) = Cs/(1.0+exp(-(x(allwn)-L1(allwn))/s)) -
         Cs/(1.0+exp(-(x(allwn)-L2(allwn))/s))+Cbkg;
       break;

    case KEY_ILS_SIGMOID_LORENTZ:
       x(allwn) = wn(allwn) - wn_center;
       Cbkg = 0.0;
       Cs = 1.0;
       if (num_para > 5)
          throw Exception("Error ! Number of ILS Parameters > 5");
       a = para(0);    // FWHM Sigmoid
       s = para(1);    // slope
       c = para(2);    // FWHM Lorentz
       d = para(3);    // Lorentz Peak 
       t = para(4);    // threshold (for GOSAT specs, .03, .006, .006)

       L1(0) = -a/2.0;
       L2(0) = a/2.0;

       Peak = Cs/(1.0+exp(L1(0)/s)) - Cs/(1.0+exp(L2(0)/s));

       ils_out(allwn) = Cs/(1.0+exp(-(x(allwn)-L1(0))/s)) -
         Cs/(1.0+exp(-(x(allwn)-L2(0))/s));
       
       for (i=0;i < wn_length;i++)
         if (ils_out(i) < t*Peak)
	   ils_out(i) = (d*c*c/4.0) / (x(i)*x(i)+(c/2.0)*(c/2.0));
       break;
       
    default:
       throw Exception("Unknown ILS function selected");
  }
}
*/

#endif
