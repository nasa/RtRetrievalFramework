// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "rt_atmosphere.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}

%base_import(observer)
%base_import(state_vector)
%import "array_ad_with_unit.i"
%import "auto_derivative.i"
%import "ground.i"

%fp_shared_ptr(FullPhysics::RtAtmosphere)
namespace FullPhysics {
  class RtAtmosphere;
}

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::RtAtmosphere>);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::RtAtmosphere>);

namespace FullPhysics {
%template(ObservableRtAtmosphere) Observable<RtAtmosphere>;
%template(ObserverRtAtmosphere) Observer<RtAtmosphere>;

class RtAtmosphere : virtual public StateVectorObserver, 
		 public Observable<RtAtmosphere> {
public:
  virtual ~RtAtmosphere();
  virtual void add_observer(Observer<RtAtmosphere>& Obs);
  virtual void remove_observer(Observer<RtAtmosphere>& Obs);
  %python_attribute(timer_info, std::string)
  %python_attribute_abstract(number_layer, int)
  %python_attribute_abstract(number_spectrometer, int)
  virtual ArrayAdWithUnit<double, 1> altitude(int spec_index) 
    const = 0;
  virtual AutoDerivative<double> 
  column_optical_depth(double wn, int spec_index, const std::string& Gas_name) const = 0;
  virtual ArrayAd<double, 1> 
    optical_depth_wrt_iv(double wn, int spec_index) const = 0;
  virtual ArrayAd<double, 1> 
    single_scattering_albedo_wrt_iv(double wn, int spec_index) const = 0;
  virtual ArrayAd<double, 3>
  scattering_moment_wrt_iv(double wn, int spec_index, int nummom = -1, 
		    int numscat = -1) const = 0;
  virtual ArrayAd<double, 1> 
  optical_depth_wrt_iv(double wn, int spec_index,
		const ArrayAd<double, 2>& iv) const = 0;
  virtual ArrayAd<double, 1> 
  single_scattering_albedo_wrt_iv(double wn, int spec_index,
			   const ArrayAd<double, 2>& iv) const = 0;
  virtual ArrayAd<double, 3>
  scattering_moment_wrt_iv(double wn, int spec_index, 
			   const ArrayAd<double, 2>& iv,
			   int nummom = -1, int numscat = -1) const = 0;
  ArrayAd<double, 1> 
    optical_depth_wrt_state_vector(double wn, int spec_index) const;
  ArrayAd<double, 1> 
    single_scattering_albedo_wrt_state_vector(double wn, int spec_index) const;
  ArrayAd<double, 3>
  scattering_moment_wrt_state_vector(double wn, int spec_index, 
     int nummom = -1, int numscat = -1) const;
  virtual ArrayAd<double, 2>
    intermediate_variable(double wn, int spec_index) const = 0;
  %python_attribute_abstract(ground, boost::shared_ptr<Ground>)
  %python_attribute_abstract(uplooking, bool)
  virtual void reset_timer();
  std::string print_to_string() const;
};
}

