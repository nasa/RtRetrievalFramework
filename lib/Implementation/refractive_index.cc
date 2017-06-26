#include "refractive_index.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

SimpleRefractiveIndex::SimpleRefractiveIndex(double bwpar, 
					     const Array<AutoDerivative<double>, 1>& press, 
					     const Array<AutoDerivative<double>, 1>& temp)
  : bwpar_(bwpar), press_(press), temp_(temp) 
{
  if(press_.rows() != temp_.rows()) {
    Exception err;
    err << "Number of pressure layers " << press_.rows() 
	<< " does not match number of temperature layers " << temp_.rows();
    throw err;
  }
}


AutoDerivative<double> SimpleRefractiveIndex::at_level(int level_index) const
{
  range_check(level_index, 0, number_level());
  return refr_index_rs(bwpar_, press_(level_index), temp_(level_index));
}

AutoDerivative<double> SimpleRefractiveIndex::at_layer(int layer_index, const AutoDerivative<double>& interp_val) const
{
  range_check(layer_index, 0, number_layer());
  range_check(value(interp_val), 0.0, 1.0);
  
  AutoDerivative<double> fp = exp( interp_val * log(press_(layer_index)) + (1 - interp_val) * log(press_(layer_index+1)) );
  AutoDerivative<double> ft = interp_val * temp_(layer_index) + (1 - interp_val) * temp_(layer_index+1);
  return refr_index_rs(bwpar_, fp, ft);
}

//-----------------------------------------------------------------------

OcoRefractiveIndex::OcoRefractiveIndex(double reference_wavelength,
				       const Array<AutoDerivative<double>, 1>& press, 
				       const Array<AutoDerivative<double>, 1>& temp,
				       const Array<AutoDerivative<double>, 1>& co2_vmr,
				       const Array<AutoDerivative<double>, 1>& h2o_vmr)
  : ref_wvl_(reference_wavelength), press_(press), temp_(temp), co2_vmr_(co2_vmr), h2o_vmr_(h2o_vmr)
{
  if( press_.rows() != temp_.rows() || 
      press_.rows() != co2_vmr_.rows() || 
      press_.rows() != h2o_vmr_.rows() ) {
    Exception err;
    err << "Number of pressure layers " << press_.rows() 
	<< " and number of temperature layers " << temp_.rows() 
	<< " and number of CO2 layers " << co2_vmr_.rows() 
	<< " and number of H2O layers " << h2o_vmr_.rows() 
	<< " must all be the same size.";
    throw err;
  }
}


AutoDerivative<double> OcoRefractiveIndex::at_level(int level_index) const
{
  range_check(level_index, 0, number_level());
  return refr_index_vn(ref_wvl_, press_(level_index), temp_(level_index), co2_vmr_(level_index), h2o_vmr_(level_index));
}

AutoDerivative<double> OcoRefractiveIndex::at_layer(int layer_index, const AutoDerivative<double>& interp_val) const
{
  range_check(layer_index, 0, number_layer());
  range_check(value(interp_val), 0.0, 1.0);

  AutoDerivative<double> fp = exp( interp_val * log(press_(layer_index)) + (1 - interp_val) * log(press_(layer_index+1)) );
  AutoDerivative<double> ft = interp_val * temp_(layer_index) + (1 - interp_val) * temp_(layer_index+1);
  AutoDerivative<double> fc = interp_val * co2_vmr_(layer_index) + (1 - interp_val) * co2_vmr_(layer_index+1);
  AutoDerivative<double> fh = interp_val * h2o_vmr_(layer_index) + (1 - interp_val) * h2o_vmr_(layer_index+1);

  return refr_index_vn(ref_wvl_, fp, ft, fc, fh);
}

bool OcoRefractiveIndex::wavelength_in_bounds(double wavelength)
{
  if ( ((wavelength >= 0.755e0) && (wavelength <= 0.785e0)) ||
       ((wavelength >= 1.58e0) && (wavelength <= 1.65e0)) || 
       ((wavelength >= 2.03e0) && (wavelength <= 2.09e0)) )
    return true;
  else
    return false;
}

//-----------------------------------------------------------------------

AutoDerivative<double> FullPhysics::refr_index_rs(double bwpar, const AutoDerivative<double>& press, const AutoDerivative<double>& temp) {
  // Original method used by ChapmanBoa from Rob Spurr
  //
  // Inputs: 
  // bwpar: base refr index param w/o 1.0 part, ie 0.000288;
  // press: pressure, Pascals
  // temp: temperature, Kelvin

  //   Standard temperature (K) and pressure (mbar).
  const double T_STANDARD = 273.16E0;
  const double P_STANDARD = 1013.25E0;
  const double STP_RATIO = T_STANDARD / P_STANDARD;

  //  Simple approximatio
  // 1.0e-2 converts pressure from Pa to mb
  double alpha = bwpar * STP_RATIO;
  AutoDerivative<double> refindex = 1.0E0 + alpha * press * 1.0e-2 / temp;

  
  return refindex;
} // End of refindex


AutoDerivative<double> FullPhysics::refr_index_vn(double wl, const AutoDerivative<double>& p, const AutoDerivative<double>& T, const AutoDerivative<double>& xc, const AutoDerivative<double>& xw) {
  // Improved refractive index method from Vijay Natraj,
  // although limited to OCO bands
  //
  //  Inputs (wl: microns, p: Pascals, T: Kelvin, xc: vmr, xw: vmr)

  //  Output
  AutoDerivative<double> nair;
  nair = 0.e0;

  AutoDerivative<double> tc = T - 273.15e0; // temperature in degrees Celsius

  if ((wl >= 0.755e0) && (wl <= 0.785e0)) {
    double wlm2 = 1.e0/  pow(wl, 2);

    //  Refractive index of standard dry air (400 ppm CO2) at 101325 Pa, 20 C
    //  Zhang et al., Appl. Opt., 47(17), 3143-3151, 2008
    //  Eq. (14)
    double nsm1 = (8015.514e0 + 2368616.e0 / (128.7459e0 - wlm2) + 
		   19085.73e0 / (50.01974e0 - wlm2)) * 1.e-8;

    //  Air density factor
    AutoDerivative<double> Dair = calc_Dair(p, tc);

    //  Refractive index of dry air (CO2 vmr = xc) at p, T
    //  Zhang et al., Appl. Opt., 47(17), 3143-3151, 2008
    //  Eq. (13), with modification for refractive index of CO2
    //  using Eq. (31)
    AutoDerivative<double> na = 1.e0 + nsm1 * (1.e0 + 0.5294e0 * (xc - 0.0004e0)) * 
      Dair / 94449.94e0;

    //  Lorentz-Lorenz factor for dry air at p, T, xc
    AutoDerivative<double> La = (pow(na,2) - 1.e0) / (pow(na,2) + 2.e0);

    //  Refractive index of pure water vapor at 1333 Pa, 20 C
    double nw = calc_nw(wlm2);

    //  Lorentz-Lorenz factor for pure water vapor at 1333 Pa, 20 C
    double Lw = (pow(nw,2) - 1.e0) / (pow(nw,2) + 2.e0);

    //  Compressibility factor for moist air (H2O vmr = xw) at p, T
    AutoDerivative<double> Z = calc_Z(p, T, tc, xw);

    //  Compressibility factor for pure water vapor at 1333 Pa, 20 C
    AutoDerivative<double> Zw = calc_Z(1333.e0, 293.15e0, 20.e0, 1.e0);

    //  Density of water vapor component of moist air at p, T, xw
    AutoDerivative<double> rhow = calc_rho(p, xw, 0.018015e0, Z, T);

    //  Density of pure water vapor at 1333 Pa, 20 C
    AutoDerivative<double> rhows = calc_rho(1333.e0, 1.e0, 0.018015e0, Zw, 293.15e0);

    //  Lorentz-Lorentz factor for moist air at p, T, xc, xw
    AutoDerivative<double> L = La + (rhow / rhows) * Lw;

    //  Refractive index of moist air at p, T, xc, xw
    nair = sqrt((1.e0 + 2.e0 * L)/(1.e0 - L));
    
  } else if (((wl >= 1.58e0) && (wl <= 1.65e0)) || 
	     ((wl >= 2.03e0) && (wl <= 2.09e0))) {

    double nw, na;
    if ((wl >= 1.58e0) && (wl <= 1.6e0)) {

      //  Refractive index of standard dry air (450 ppm CO2) at 101325 Pa, 293.16 K
      //  Colavita et al., Publ. Astron. Soc. Pac., 116, 876-885, 2004
      //  Table 1

      double n1 = 2.6859e-4;
      double n2 = 2.6855e-4;
      na = 1.e0 + n1 + ((wl - 1.55e0) / 0.05e0) * (n2 - n1);

      //  Refractive index of pure water vapor at 410 Pa, 293.16 K
      //  Colavita et al., Publ. Astron. Soc. Pac., 116, 876-885, 2004
      //  Table 1
      n1 = 9.164e-7;
      n2 = 9.154e-7;
      nw = 1.e0 + n1 + ((wl - 1.55e0) / 0.05e0) * (n2 - n1);

    } else if ((wl > 1.6e0) && (wl <= 1.65e0)) {
      double n1 = 2.6855e-4;
      double n2 = 2.6851e-4;
      na = 1.e0 + n1 + ((wl - 1.6e0) / 0.05e0) * (n2 - n1);
      
      n1 = 9.154e-7;
      n2 = 9.144e-7;
      nw = 1.e0 + n1 + ((wl - 1.6e0) / 0.05e0) * (n2 - n1);

    } else if ((wl >= 2.03e0) && (wl <= 2.05e0)) {
      double n1 = 2.6833e-4;
      double n2 = 2.6831e-4;
      na = 1.e0 + n1 + ((wl - 2.e0) / 0.05e0) * (n2 - n1);

      n1 = 9.092e-7;
      n2 = 9.076e-7;
      nw = 1.e0 + n1 + ((wl - 2.e0) / 0.05e0) * (n2 - n1);

    } else if ((wl > 2.05e0) && (wl <= 2.09e0)) {
      double n1 = 2.6831e-4;
      double n2 = 2.683e-4;
      na = 1.e0 + n1 + ((wl - 2.05e0) / 0.05e0) * (n2 - n1);

      n1 = 9.076e-7;
      n2 = 9.061e-7;
      nw = 1.e0 + n1 + ((wl - 2.05e0) / 0.05e0) * (n2 - n1);

    }

    //  Compressibility factor for moist air at p, T, xw
    AutoDerivative<double> Z = calc_Z(p, T, tc, xw);

    //  Compressibility factor for standard dry air at 101325 Pa, 293.16 K
    AutoDerivative<double> Za = calc_Z(101325.e0, 293.16e0, 20.01e0, 0.e0);

    //  Compressibility factor for pure water vapor at 410 Pa, 20.01 C
    AutoDerivative<double> Zw = calc_Z(410.e0, 293.16e0, 20.01e0, 1.e0);

    //  Molar mass of dry air, xc
    AutoDerivative<double> Ma = (28.9635e0 + 12.011e0 * (xc - 0.0004e0)) * 1.e-3;

    //  Density of dry component of moist air at p, T, xc, xw
    AutoDerivative<double> rhoa = calc_rho(p, 1.e0-xw, Ma, Z, T);

    //  Molar mass of standard dry air, 450 ppm CO2
    double Mas = (28.9635e0 + 12.011e0 * (0.00045e0 - 0.0004e0)) * 1.e-3;

    //  Density of standard dry air at 101325 Pa, 293.16 K
    AutoDerivative<double> rhoas = calc_rho(101325.e0, 1.e0, Mas, Za, 293.16e0);

    //  Density of water vapor component of moist air at p, T, xw
    AutoDerivative<double> rhow = calc_rho(p, xw, 0.018015e0, Z, T);
    
    //  Density of pure water vapor at 1333 Pa, 20 C
    AutoDerivative<double> rhows = calc_rho(1333.e0, 1.e0, 0.018015e0, Zw, 293.15e0);
    
    //  Refractive index of moist air at p, T, xc, xw
    nair = 1.e0 + (rhoa / rhoas) * (na - 1.e0) + (rhow / rhows) * (nw - 1.e0);
      
  } else {    
    throw Exception("Can not compute refractive index outside spectral range handled by this function");
      
  }

  return nair;

}

AutoDerivative<double> FullPhysics::calc_Dair(const AutoDerivative<double>& p, const AutoDerivative<double>& tc) {

  //  Zhang et al., Appl. Opt., 47(17), 3143-3151, 2008
  //  Eq. (11)
  AutoDerivative<double> Dair =
    p * (1.e0 + p * (0.621811e0 - 0.0126531e0 * tc + 0.000066e0 * pow(tc,2)) *
	 1.e-8) / (1.e0 + 0.003661e0 * tc);

  return Dair;

}

double FullPhysics::calc_nw(double wlm2) {

  //  Local variables
  double w0, w1, w2, w3;

  //  Ciddor, Appl. Opt., 35(9), 1566-1573, 1996
  //  Eq. (3)
  w0 = 295.235e0;
  w1 = 2.6422;
  w2 = -0.03238e0;
  w3 = 0.004028e0;
  double nw = 1.022e0*(w0+w1*wlm2+w2*pow(wlm2,2)+w3*pow(wlm2,3))*1.e-8+1.e0;

  return nw;

}

AutoDerivative<double> FullPhysics::calc_Z(const AutoDerivative<double>& p, const AutoDerivative<double>& T, const AutoDerivative<double>& tc, const AutoDerivative<double>& xw) {

  //  Local variables
  double a0, a1, a2, b0, b1, c0, c1, d, e;

  //  Ciddor, Appl. Opt., 35(9), 1566-1573, 1996
  //  Eq. (12)
  a0 = 1.58123e-6;
  a1 = -2.9331e-8;
  a2 = 1.1043e-10;
  b0 = 5.707e-6;
  b1 = -2.051e-8;
  c0 = 1.9898e-4;
  c1 = -2.376e-6;
  d = 1.83e-11;
  e = -0.765e-8;
  AutoDerivative<double> Z = 1.e0-(p/T)*(a0+a1*tc+a2*pow(tc,2)+(b0+b1*tc)*xw+(c0+c1*tc)*pow(xw,2))+ 
    pow((p/T),2)*(d+e*pow(xw,2));

  return Z;
}


AutoDerivative<double> FullPhysics::calc_rho(const AutoDerivative<double>& p, const AutoDerivative<double>& x, const AutoDerivative<double>& M, const AutoDerivative<double>& Z, const AutoDerivative<double>& T) {

  const double R = 8.314510e0;

  //  Ciddor, Appl. Opt., 41(12), 2292-2298, 2002
  //  Eq. (3)
  AutoDerivative<double> rho = p*x*M/(Z*R*T);

  return rho;

}
