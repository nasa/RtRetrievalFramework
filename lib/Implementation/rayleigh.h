#ifndef RAYLEIGH_H
#define RAYLEIGH_H
#include "pressure.h"
#include "altitude.h"
#include "constant.h"
#include "default_constant.h"
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This class calculates the Rayleigh portion of the optical depth
*******************************************************************/
class Rayleigh: public Observer<Pressure>, public Observer<Altitude>,
		  public Printable<Rayleigh> {
public:
  Rayleigh(const boost::shared_ptr<Pressure>& Pres, 
	   const std::vector<boost::shared_ptr<Altitude> >& Alt,
	   const Constant& C);
  
  virtual void notify_update(const Pressure& P)
  { cache_is_stale = true;}
  virtual void notify_update(const Altitude& A)
  {cache_is_stale = true;}

  ArrayAd<double, 1> optical_depth_each_layer(double wn, int spec_index) const;
  static DoubleWithUnit cross_section(const DoubleWithUnit& W,
				      const Constant& C = DefaultConstant());
  virtual void print(std::ostream& Os) const
  { Os << "Rayleigh"; }
private:
  boost::shared_ptr<Pressure> pres;
  std::vector<boost::shared_ptr<Altitude> > alt;
  mutable bool cache_is_stale;
  mutable ArrayAd<double, 2> part_independent_wn;
  void fill_cache() const;
  // Constants. We get this from the Constant class, but stash a copy
  // of them here.
  double a, b, depolar_fact, molar_weight_dry_air;
};
}
#endif
