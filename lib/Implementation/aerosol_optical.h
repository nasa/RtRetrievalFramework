#ifndef AEROSOL_OPTICAL_H
#define AEROSOL_OPTICAL_H
#include "aerosol.h"
#include "aerosol_property.h"
#include "aerosol_extinction.h"
#include <limits>
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  Implementation of Aerosol. This particular implementation does 
  the aerosol calculation by using the aerosol optical properties.
*******************************************************************/
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
  virtual ~AerosolOptical() {}
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
  virtual void notify_update(const StateVector& Sv) 
  {
    nvar = Sv.state_with_derivative().number_variable();
    notify_update_do(*this); 
  }

  virtual ArrayAd<double, 2> optical_depth_each_layer(double wn) 
    const;
  virtual ArrayAd<double, 1> 
  ssa_each_layer(double wn, int particle_index,
		 const ArrayAd<double, 1>& Od) const;
  virtual ArrayAd<double, 1> 
  ssa_each_layer(double wn) const;

//-----------------------------------------------------------------------
/// For performance, we cache some data as we calculate it. This
/// becomes stale when the pressure is changed, so we observe press
/// and when it changes notify other observers that we have changed.
//-----------------------------------------------------------------------

  virtual void notify_update(const Pressure& P)
  {
    cache_is_stale = true;
    notify_update_do(*this);
  }
  virtual void notify_update(const AerosolExtinction& A)
  {
    cache_is_stale = true;
    notify_update_do(*this);
  }
  virtual void notify_update(const AerosolProperty& A)
  {
    cache_is_stale = true;
    notify_update_do(*this);
  }

  virtual ArrayAd<double, 3> pf_mom(double wn, int pindex) const;
  virtual blitz::Array<double, 3> pf_mom(double wn, 
			    const blitz::Array<double, 2>& frac_aer) const;
  virtual ArrayAd<double, 3> pf_mom(double wn, 
				    const ArrayAd<double, 2>& frac_aer,
				    int nummom = -1, int numscat = -1) const;
  virtual int number_particle() const { return (int) aext.size(); }

  virtual void print(std::ostream& Os) const;

//-----------------------------------------------------------------------
// Optical depth of each aerosol type or total
//-----------------------------------------------------------------------
  double aerosol_optical_depth
  (int aer_idx,
   double pmin = std::numeric_limits<double>::min(),
   double pmax = std::numeric_limits<double>::max()) const;
  double aerosol_optical_depth_total
  (double pmin = std::numeric_limits<double>::min(),
   double pmax = std::numeric_limits<double>::max()) const;

//-----------------------------------------------------------------------
/// Clone a Aerosol object. Note that the cloned version will *not*
/// be attached to and StateVector or Observer<Pressure>, although you
/// can of course attach them after receiving the cloned object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Aerosol> clone() const 
  { return clone(press->clone(), rh->clone()); }
  virtual boost::shared_ptr<Aerosol> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<RelativeHumidity>& Rh) const;
  std::vector<std::string> aerosol_name() const;

  blitz::Array<std::string, 1> aerosol_name_arr() const;

//-----------------------------------------------------------------------
/// Return aerosol extinction.
//-----------------------------------------------------------------------
  const boost::shared_ptr<AerosolExtinction>& aerosol_extinction(int i) const
  { range_check(i, 0, number_particle()); return aext[i];}

//-----------------------------------------------------------------------
/// Set AerosolExtinction.
//-----------------------------------------------------------------------

  void aerosol_extinction(int i, const boost::shared_ptr<AerosolExtinction>& V)
  { range_check(i, 0, number_particle()); aext[i] = V; 
    cache_is_stale = true; notify_update_do(*this); }

//-----------------------------------------------------------------------
/// Return aerosol property
//-----------------------------------------------------------------------

  const boost::shared_ptr<AerosolProperty>& aerosol_property(int i) const
  { range_check(i, 0, number_particle()); return aprop[i];}

//-----------------------------------------------------------------------
/// Set AerosolProperty.
//-----------------------------------------------------------------------

  void aerosol_property(int i, const boost::shared_ptr<AerosolProperty>& V)
  { range_check(i, 0, number_particle()); aprop[i] = V; 
    cache_is_stale = true; notify_update_do(*this); }

//-----------------------------------------------------------------------
/// Return pressure
//-----------------------------------------------------------------------

  const boost::shared_ptr<Pressure>& pressure() const
  { return press; }
private:
  std::vector<boost::shared_ptr<AerosolExtinction> > aext;
  std::vector<boost::shared_ptr<AerosolProperty> > aprop;
  boost::shared_ptr<Pressure> press;
  boost::shared_ptr<RelativeHumidity> rh;
  double reference_wn_;
  // We cache to part of optical_depth_each_layer calculation that
  // independent of wn.
  mutable bool cache_is_stale;
  mutable ArrayAd<double, 2> od_ind_wn;
  mutable int nlay;
  void fill_cache() const;
  int nvar;
};
}
#endif
