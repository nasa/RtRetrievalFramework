#ifndef HRES_WRAPPER_H
#define HRES_WRAPPER_H
#include "radiative_transfer_single_wn.h"

namespace FullPhysics {
/****************************************************************//**
  For timing purposes, it can be useful to separate out the high
  resolution radiative transfer vs. the low resolution. Valgrind
  lumps them all together, since it is organized by function call
  and the same function is called for both high resolution and low
  resolution.

  This class provides an easy work around. It just forwards 
  everything to an actual RadiativeTransferSingleWn class, but it
  wraps this in a separate function call. This then means valgrind
  lists this as separate.
*******************************************************************/

class HresWrapper : public RadiativeTransferSingleWn
{
public:
  HresWrapper(const boost::shared_ptr<RadiativeTransferSingleWn>& Rt)
    : RadiativeTransferSingleWn(Rt->stokes_coefficient(), 
				Rt->atmosphere_ptr()),
      rt_(Rt) {}
  virtual ~HresWrapper() {}
  virtual int number_stokes() const { return rt_->number_stokes(); }
  virtual int number_stream() const { return rt_->number_stream(); }
  virtual blitz::Array<double, 1> stokes_single_wn
  (double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn
  (double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;
  virtual boost::shared_ptr<RadiativeTransfer> rt() const {return rt_;}
private:
  boost::shared_ptr<RadiativeTransferSingleWn> rt_;
};
}
#endif
