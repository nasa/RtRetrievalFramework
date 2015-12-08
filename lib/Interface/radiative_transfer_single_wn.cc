#include "radiative_transfer_single_wn.h"
#include "ostream_pad.h"
using namespace FullPhysics;
using namespace blitz;

Array<double, 2> 
RadiativeTransferSingleWn::stokes(const SpectralDomain& Spec_domain,
				  int Spec_index) const
{
  Array<double, 1> wn(Spec_domain.wavenumber());
  boost::shared_ptr<boost::progress_display> disp = progress_display(wn);
  Array<double, 2> res(wn.rows(), number_stokes());
  for(int i = 0; i < wn.rows(); ++i) {
    res(i, Range::all()) = stokes_single_wn(wn(i), Spec_index);
    if(disp)	
      *disp += 1;
  }
  return res;
}

ArrayAd<double, 2> 
RadiativeTransferSingleWn::stokes_and_jacobian(const SpectralDomain& Spec_domain,
					       int Spec_index) const
{
  Array<double, 1> wn(Spec_domain.wavenumber());
  if(wn.rows() < 1)		// Handle degenerate case.
    return ArrayAd<double, 2>(0,number_stokes(),0);
  boost::shared_ptr<boost::progress_display> disp = progress_display(wn);
  ArrayAd<double, 1> t =
    stokes_and_jacobian_single_wn(wn(0), Spec_index);
  ArrayAd<double, 2> res(wn.rows(), number_stokes(), t.number_variable());
  res(0, Range::all()) = t;
  for(int i = 1; i < wn.rows(); ++i) {
    res(i, Range::all()) = stokes_and_jacobian_single_wn(wn(i), Spec_index);
    if(disp)	
      *disp += 1;
  }
  return res;
}

void RadiativeTransferSingleWn::print(std::ostream& Os, bool Short_form) const
{
  RadiativeTransferFixedStokesCoefficient::print(Os, Short_form);
  OstreamPad opad(Os, "  ");
  if(!Short_form) {
    Os << "\nAtmosphere:\n";
    opad << *atm << "\n";
    opad.strict_sync();
  }
}
