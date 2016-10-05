// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "ils_convolution.h"
#include "sub_state_vector_array.h"
%}
%base_import(ils)
%base_import(observer)
%base_import(dispersion)
%import "double_with_unit.i"
%import "state_vector.i"
%import "array_ad.i"
%import "ils_function.i"
%fp_shared_ptr(FullPhysics::IlsConvolution);

namespace FullPhysics {
class IlsConvolution : public Ils, public Observer<Dispersion> {
public:
  IlsConvolution(const boost::shared_ptr<Dispersion>& Disp,
		 const boost::shared_ptr<IlsFunction>& Ils_func,
		 const DoubleWithUnit& Ils_half_width = DoubleWithUnit(20, units::inv_cm));
  IlsConvolution(const boost::shared_ptr<Dispersion>& Disp,
		 const boost::shared_ptr<IlsFunction>& Ils_func,
		 double Ils_half_width);
  virtual void notify_update(const StateVector& Sv);
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
  virtual void notify_update(const Dispersion& D);
  virtual blitz::Array<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const blitz::Array<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const;
  virtual ArrayAd<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const ArrayAd<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const;
  virtual boost::shared_ptr<Ils> clone() const;
  %python_attribute(ils_function, boost::shared_ptr<IlsFunction>);
  %python_attribute(dispersion, boost::shared_ptr<Dispersion>);
};
}
