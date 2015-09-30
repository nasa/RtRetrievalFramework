// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "aerosol.h"
%}
%base_import(state_vector)
%base_import(pressure)
%base_import(aerosol_extinction)
%base_import(aerosol_property)
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
               public Observer<Pressure>,
	       public Observer<AerosolExtinction>,
	       public Observer<AerosolProperty>,
	       public Observable<Aerosol> {
public:
  std::string print_to_string() const;
  Aerosol(const std::vector<boost::shared_ptr<AerosolExtinction> >& Aext,
	  const std::vector<boost::shared_ptr<AerosolProperty> >& Aerosol_prop,
	  const boost::shared_ptr<Pressure>& Press,
	  double Reference_wn = 1e4/0.755);
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
  virtual void notify_update(const StateVector& Sv);
  virtual void add_observer(Observer<Aerosol>& Obs);
  virtual void remove_observer(Observer<Aerosol>& Obs);

  ArrayAd<double, 2> optical_depth_each_layer(double wn) 
    const;
  ArrayAd<double, 1> 
  ssa_each_layer(double wn, int particle_index,
		 const ArrayAd<double, 1>& Od) const;
  ArrayAd<double, 1> 
  ssa_each_layer(double wn) const;
  virtual void notify_update(const Pressure& P);
  virtual void notify_update(const AerosolExtinction& A);
  virtual void notify_update(const AerosolProperty& A);
  blitz::Array<double, 2> pf_mom(double wn, int pindex) const;
  blitz::Array<double, 3> pf_mom(double wn, 
				 const blitz::Array<double, 2>& frac_aer) const;
  ArrayAd<double, 3> pf_mom(double wn, 
	 const ArrayAd<double, 2>& frac_aer,
	 int nummom = -1, int numscat = -1) const;
  %python_attribute(number_particle, int)
  double aerosol_optical_depth
  (int aer_idx,
   double pmin = std::numeric_limits<double>::min(),
   double pmax = std::numeric_limits<double>::max()) const;
  double aerosol_optical_depth_total
  (double pmin = std::numeric_limits<double>::min(),
   double pmax = std::numeric_limits<double>::max()) const;
  boost::shared_ptr<Aerosol> clone() const;
  boost::shared_ptr<Aerosol> 
  clone(const boost::shared_ptr<Pressure>& Press) const;
  %python_attribute(aerosol_name, std::vector<std::string>)
  %python_attribute(aerosol_name_arr, blitz::Array<std::string, 1>)

  boost::shared_ptr<AerosolExtinction> aerosol_extinction(int i) const;
  boost::shared_ptr<AerosolProperty> aerosol_property(int i) const;
};
}
