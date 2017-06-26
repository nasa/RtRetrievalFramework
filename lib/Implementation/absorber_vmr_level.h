#ifndef ABSORBER_VMR_LEVEL_H
#define ABSORBER_VMR_LEVEL_H
#include "absorber_vmr_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This particular implementation uses the state vector values as
  the VMR for each pressure level.
*******************************************************************/
class AbsorberVmrLevel : public AbsorberVmrImpBase {
public:
  AbsorberVmrLevel(const boost::shared_ptr<Pressure>& Press,
		   const blitz::Array<double, 1>& Vmr, 
		   const blitz::Array<bool, 1>& Vmr_flag,
		   const std::string& Gas_name);
  virtual ~AbsorberVmrLevel() {}
  virtual void print(std::ostream& Os) const;
  virtual std::string state_vector_name_i(int i) const
  { return gas_name() + " VMR for Press Lvl " + 
      boost::lexical_cast<std::string>(i + 1); }
  virtual boost::shared_ptr<AbsorberVmr> clone() const
  { return clone(boost::shared_ptr<Pressure>()); }
  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const;

//-----------------------------------------------------------------------
/// VMR on the pressure grid. This is just coeff.value, but this is
/// useful for generating output.
//-----------------------------------------------------------------------
  blitz::Array<double, 1> vmr_profile() const 
  { return coeff.value(); }

//-----------------------------------------------------------------------
/// Covariance of vmr profile
//-----------------------------------------------------------------------

  blitz::Array<double, 2> vmr_covariance() const
  {
    using namespace blitz;
    firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
    ArrayAd<double, 1> vmrv(coeff);
    Array<double, 2> res(vmrv.rows(), vmrv.rows());
    if(vmrv.is_constant())
      res = 0;			// Normal only encounter this in
  				// testing, when we haven't yet set up
  				// a covariance matrix and state vector.
    else {
      
      Array<double, 2> dvmr_dstate(vmrv.jacobian());
      Array<double, 2> t(dvmr_dstate.rows(), sv_cov_full.cols()) ;
      t = sum(dvmr_dstate(i1, i3) * sv_cov_full(i3, i2), i3);
      res = sum(t(i1, i3) * dvmr_dstate(i2, i3), i3);
    }
    return res;
  }

//-----------------------------------------------------------------------
/// Uncertainty of VMR
//-----------------------------------------------------------------------

  blitz::Array<double, 1> vmr_uncertainty() const
  {
    blitz::Array<double, 2> vmrcov(vmr_covariance());
    blitz::Array<double, 1> res(vmrcov.rows());
    for(int i = 0; i < res.rows(); ++i)
      res(i) = (vmrcov(i, i) > 0 ? sqrt(vmrcov(i, i)) : 0.0);
    return res;
  }
protected:
  virtual void calc_vmr() const;
};
}
#endif
