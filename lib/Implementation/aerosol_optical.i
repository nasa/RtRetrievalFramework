// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "sub_state_vector_array.h"
#include "aerosol_optical.h"
#include "absorber.h"
#include "temperature.h"
#include "altitude.h"
%}
%base_import(state_vector)
%base_import(pressure)
%base_import(aerosol_extinction)
%base_import(aerosol_property)
%base_import(aerosol)
%import "relative_humidity.i"
%fp_shared_ptr(FullPhysics::AerosolOptical);

namespace FullPhysics {
class AerosolOptical: public Aerosol,
               public Observer<Pressure>,
	       public Observer<AerosolExtinction>,
	       public Observer<AerosolProperty> {
public:
  AerosolOptical(const std::vector<boost::shared_ptr<AerosolExtinction> >& Aext,
	  const std::vector<boost::shared_ptr<AerosolProperty> >& Aerosol_prop,
	  const boost::shared_ptr<Pressure>& Press,
	  const boost::shared_ptr<RelativeHumidity>& Rh,
	  double Reference_wn = 1e4/0.755);
  virtual ArrayAd<double, 2> optical_depth_each_layer(double wn) 
    const;
  virtual ArrayAd<double, 1> 
  ssa_each_layer(double wn, int particle_index,
		 const ArrayAd<double, 1>& Od) const;
  virtual ArrayAd<double, 1> 
  ssa_each_layer(double wn) const;
  virtual void notify_update(const Pressure& P);
  virtual void notify_update(const AerosolExtinction& A);
  virtual void notify_update(const AerosolProperty& A);
  virtual ArrayAd<double, 3> pf_mom(double wn, int pindex) const;
  virtual blitz::Array<double, 3> pf_mom(double wn, 
				 const blitz::Array<double, 2>& frac_aer) const;
  virtual ArrayAd<double, 3> pf_mom(double wn, 
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
  virtual boost::shared_ptr<Aerosol> clone() const;
  virtual boost::shared_ptr<Aerosol> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<RelativeHumidity>& Rh) const;
  %python_attribute(aerosol_name, std::vector<std::string>);
  %python_attribute(aerosol_name_arr, blitz::Array<std::string, 1>);
  %python_attribute(pressure, boost::shared_ptr<Pressure>);
  boost::shared_ptr<AerosolExtinction> aerosol_extinction(int i) const;
  void aerosol_extinction(int i, const boost::shared_ptr<AerosolExtinction>& V);
  boost::shared_ptr<AerosolProperty> aerosol_property(int i) const;
  void aerosol_property(int i, const boost::shared_ptr<AerosolProperty>& V);
};
}
