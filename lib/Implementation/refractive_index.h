#ifndef REFR_INDEX_H
#define REFR_INDEX_H

#include <blitz/array.h>
#include "auto_derivative.h"

namespace FullPhysics {

  AutoDerivative<double> refr_index_rs(double bwpar, const AutoDerivative<double>& press, const AutoDerivative<double>& temp);
  AutoDerivative<double> refr_index_vn(double wl, const AutoDerivative<double>& p, const AutoDerivative<double>& T, const AutoDerivative<double>& xc, const AutoDerivative<double>& xw);

  AutoDerivative<double> calc_Dair(const AutoDerivative<double>& p, const AutoDerivative<double>& tc);
  double calc_nw(double wlm2);
  AutoDerivative<double> calc_Z(const AutoDerivative<double>& p, const AutoDerivative<double>& T, const AutoDerivative<double>& tc, const AutoDerivative<double>& xw);
  AutoDerivative<double> calc_rho(const AutoDerivative<double>& p, const AutoDerivative<double>& x, const AutoDerivative<double>& M, const AutoDerivative<double>& Z, const AutoDerivative<double>& T);

/****************************************************************//**
  This class provides an interface for classes that compute a
  refractive index for an atmosphere at a given level or layer.
*******************************************************************/

class AtmRefractiveIndex { 

public:

  virtual ~AtmRefractiveIndex() { }
  //-----------------------------------------------------------------------
  /// Number of levels of atmospheric data represented
  //-----------------------------------------------------------------------
  virtual int number_level() const = 0;

  //-----------------------------------------------------------------------
  /// Number of layers of atmospheric data represented
  //-----------------------------------------------------------------------
  virtual int number_layer() const { return number_level() - 1;  }

  //-----------------------------------------------------------------------
  /// Returns refractive index at a level index
  //-----------------------------------------------------------------------
  virtual AutoDerivative<double> at_level(int level_index) const = 0;

  //-----------------------------------------------------------------------
  /// Returns refractive index at all levels
  //-----------------------------------------------------------------------
  virtual blitz::Array<AutoDerivative<double>, 1> level_values() const {
    blitz::Array<AutoDerivative<double>, 1> res(number_level());
    for(int li = 0; li < number_level(); li++)
      res(li) = at_level(li);
    return res;
  }

  //-----------------------------------------------------------------------
  /// Returns refractive index at a layer index interpolated between
  /// bounding levels using a linear interpolation factor [0,1]
  //-----------------------------------------------------------------------
  virtual AutoDerivative<double> at_layer(int layer_index, const AutoDerivative<double>& interp_val) const = 0;

  //-----------------------------------------------------------------------
  /// Returns refractive index at all layers according to a mapping
  /// of where each layer value should be calculated between levels
  //-----------------------------------------------------------------------
  virtual blitz::Array<AutoDerivative<double>, 1> layer_values(const blitz::Array<AutoDerivative<double>, 1>& interp_vals) const {
    range_max_check(interp_vals.rows(), number_layer()+1);
    blitz::Array<AutoDerivative<double>, 1> res(number_layer());
    for(int li = 0; li < number_layer(); li++)
      res(li) = at_layer(li, interp_vals(li));
    return res;
  }

  //-----------------------------------------------------------------------
  /// Convenience routine to returns refractive index ay layer index 
  /// where layer data is taken as mean of bounding level values.
  //-----------------------------------------------------------------------
  virtual AutoDerivative<double> at_layer_midpoint(int layer_index) const {
    return at_layer(layer_index, 0.5);
  }

  //-----------------------------------------------------------------------
  /// Returns refractive index at the middle of each layer
  //-----------------------------------------------------------------------
  virtual blitz::Array<AutoDerivative<double>, 1> layer_midpoint_values() const {
    blitz::Array<AutoDerivative<double>, 1> interp_vals(number_layer());
    interp_vals = 0.5;
    return layer_values(interp_vals);
  }

};

/****************************************************************//**
  Computes refractive index per level/layer using a simple
  approximation using pressure and temperature.
*******************************************************************/

class SimpleRefractiveIndex : public AtmRefractiveIndex { 

public:  
  virtual ~SimpleRefractiveIndex() {}
  SimpleRefractiveIndex(double bwpar, 
			const blitz::Array<AutoDerivative<double>, 1>& press, 
			const blitz::Array<AutoDerivative<double>, 1>& temp);
  
  virtual int number_level() const { return press_.rows(); }
  virtual AutoDerivative<double> at_level(int level_index) const;
  virtual AutoDerivative<double> at_layer(int layer_index, const AutoDerivative<double>& interp_val) const;

private:
  double bwpar_;
  blitz::Array<AutoDerivative<double>, 1> press_;
  blitz::Array<AutoDerivative<double>, 1> temp_;

};

/****************************************************************//**
  Computes refractive index per level/layer using a method that
  takes into account knowledge of OCO spectral domains.

  However, this method will fail when using spectral ranges
  outside of the OCO bands.
*******************************************************************/

class OcoRefractiveIndex : public AtmRefractiveIndex { 

public:  

  /// Creates a refractive index class that can better calculate
  /// refractive index for data with in the OCO spectral bounds
  /// reference_wavelength - microns
  /// press - Pascals
  /// temp - Kelvin
  /// co2_vmr, h2o_vmr - VMR
  OcoRefractiveIndex(double reference_wavelength, 
		     const blitz::Array<AutoDerivative<double>, 1>& press, 
		     const blitz::Array<AutoDerivative<double>, 1>& temp,
		     const blitz::Array<AutoDerivative<double>, 1>& co2_vmr,
		     const blitz::Array<AutoDerivative<double>, 1>& h2o_vmr);
  virtual ~OcoRefractiveIndex() {}
  virtual int number_level() const { return press_.rows(); }
  virtual AutoDerivative<double> at_level(int level_index) const;
  virtual AutoDerivative<double> at_layer(int layer_index, const AutoDerivative<double>& interp_val) const;

  /// Determines if the passed wavelength value in microns is within the
  /// bounds of those accepted by the class
  static bool wavelength_in_bounds(double wavelength);

private:
  double ref_wvl_;
  blitz::Array<AutoDerivative<double>, 1> press_;
  blitz::Array<AutoDerivative<double>, 1> temp_;
  blitz::Array<AutoDerivative<double>, 1> co2_vmr_;
  blitz::Array<AutoDerivative<double>, 1> h2o_vmr_;
};


}

#endif
