// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "sub_state_vector_array.h"
#include "aerosol.h"
#include "absorber.h"
#include "temperature.h"
#include "altitude.h"
%}
%base_import(state_vector)
%base_import(pressure)
%import "sub_state_vector_array.i"
%import "absorber.i"
%import "relative_humidity.i"
%fp_shared_ptr(FullPhysics::Aerosol);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Aerosol>)
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Aerosol>)

namespace FullPhysics {
class Aerosol;
}

namespace FullPhysics {
%template(ObserverAerosol) FullPhysics::Observer<Aerosol>;
%template(ObservableAerosol) FullPhysics::Observable<Aerosol>;
class Aerosol: public StateVectorObserver,
	       public Observable<Aerosol> {
public:
  std::string print_to_string() const;
  virtual void add_observer(Observer<Aerosol>& Obs);
  virtual void remove_observer(Observer<Aerosol>& Obs);

  virtual ArrayAd<double, 3> pf_mom(double wn, 
         const ArrayAd<double, 2>& frac_aer,
         int nummom = -1, int numscat = -1) const = 0;
  virtual ArrayAd<double, 2> optical_depth_each_layer(double wn) 
    const = 0;
  virtual ArrayAd<double, 1> 
  ssa_each_layer(double wn, int particle_index,
		 const ArrayAd<double, 1>& Od) const = 0;

  %python_attribute(number_particle, int)
  boost::shared_ptr<Aerosol> clone() const;
  boost::shared_ptr<Aerosol> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<RelativeHumidity>& Rh) const;
};
}
