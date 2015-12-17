// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "rayleigh.h"
#include "temperature.h"
%}
%base_import(pressure)
%base_import(temperature)
%base_import(altitude)
%import "constant.i"
%import "double_with_unit.i"
%fp_shared_ptr(FullPhysics::Rayleigh);
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Rayleigh>);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Rayleigh>);

namespace FullPhysics {
  class Rayleigh;
%template(ObservableRayleigh) FullPhysics::Observable<Rayleigh>;
%template(ObserverRayleigh) FullPhysics::Observer<Rayleigh>;

class Rayleigh: public Observer<Pressure>, public Observer<Altitude> {
public:
  Rayleigh(const boost::shared_ptr<Pressure>& Pres, 
	   const std::vector<boost::shared_ptr<Altitude> >& Alt,
   	   const Constant& C);
  virtual void notify_update(const Pressure& P);
  virtual void notify_update(const Altitude& A);
  ArrayAd<double, 1> optical_depth_each_layer(double wn,
					      int spec_index) const;
  static DoubleWithUnit cross_section(const DoubleWithUnit& W);
  static DoubleWithUnit cross_section(const DoubleWithUnit& W,
				      const Constant& C);
};
}
