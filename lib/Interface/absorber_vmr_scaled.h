#ifndef ABSORBER_VMR_SCALED_H
#define ABSORBER_VMR_SCALED_H
#include "absorber_vmr_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This interface encapsulates the usage of a vmr profile along
  along with a scale factor which is retrieved.
*******************************************************************/
class AbsorberVmrScaled : public AbsorberVmrImpBase {
public:
  AbsorberVmrScaled(const boost::shared_ptr<Pressure>& Press,
		    double Scale,                         
		    bool Scale_flag,
		    const std::string& Gas_name);

  virtual ~AbsorberVmrScaled() {}
  virtual void print(std::ostream& Os) const;

  virtual std::string state_vector_name_i(int i) const
  { return gas_name() + " Scaling factor"; }

  virtual boost::shared_ptr<AbsorberVmr> clone() const
  { return clone(press->clone()); }

  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const = 0;

  //-----------------------------------------------------------------------
  /// Scale factor.
  //-----------------------------------------------------------------------

  double scale_factor() const { return coeff(0).value(); }

  //-----------------------------------------------------------------------
  /// Uncertainty of scale factor.
  //-----------------------------------------------------------------------

  double scale_uncertainty() const
  { return (sv_cov_sub.rows() > 0 && sv_cov_sub(0,0) > 0 ? 
	    sqrt(sv_cov_sub(0,0)) : 0); }

  //-----------------------------------------------------------------------
  /// VMR profile used with the offset
  //-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> vmr_profile() const = 0;

  //-----------------------------------------------------------------------
  /// Pressure levels that absorber vmr is on
  //-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> pressure_profile() const = 0;

protected:
  virtual void calc_vmr() const;
};
}
#endif
