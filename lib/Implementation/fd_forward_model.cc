#include "fd_forward_model.h"
using namespace FullPhysics;
using namespace blitz;

#include "logger.h"

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

FdForwardModel::FdForwardModel(const boost::shared_ptr<ForwardModel>& Real_forward_model,
			       const boost::shared_ptr<StateVector>& Sv,
			       const blitz::Array<double, 1>& Perturbation)
  : real_fm(Real_forward_model), statev(Sv), perturb(Perturbation)
{
}

// See base class for description.

Spectrum FdForwardModel::radiance(int Spec_index, bool Skip_jacobian) const
{
  if(Skip_jacobian)
    return real_fm->radiance(Spec_index, Skip_jacobian);

  // Save initial state of statevector to set back at end of jacobian looping
  Array<double, 1> initial_sv( statev->state() );

  // Retrieve unperturbed radiance value
  Logger::info() << "Finite Difference FM: Unperturbed radiances\n";
  Spectrum rad(real_fm->radiance(Spec_index, true));

  // Set up result ArrayAD class
  ArrayAd<double, 1> res(rad.spectral_range().data().rows(), 
			 statev->state_with_derivative().number_variable());
  res.value() = rad.spectral_range().data();

  // Make sure statevector size is same as perturbation array
  if( initial_sv.extent(firstDim) != perturb.extent(firstDim) ) {
    Exception err;
    err << "Statevector size does not match size of perturbation array, can not do FD jacobian calculations";
    throw(err);
  }

  // Loop over statevector perturbing each item in turn, saving value into jacobian array
  Array<double, 1> current_sv(initial_sv.extent(firstDim));
  for( int sv_idx = 0; sv_idx < initial_sv.extent(firstDim); sv_idx++) {
    Logger::info() << "Finite Difference FM: Perturbation of state vector element #" << sv_idx+1 << "\n";
    if( perturb(sv_idx) > 0.0 ) {
      // Copy original and add perturbation for current element
      current_sv = initial_sv;
      current_sv(sv_idx) += perturb(sv_idx);
    
      // Update statevector so calculations are done w/ perturbed values
      statev->update_state(current_sv);

      // Calculate FD jacobian value      
      res.jacobian()(Range::all(), sv_idx) = 
	(real_fm->radiance(Spec_index, true).spectral_range().data() - 
	 rad.spectral_range().data()) / perturb(sv_idx);;

    } else {
      // Set to all zeros if perturbation is 0.0 since 
      // the perturbation would have no effect
      res.jacobian()(Range::all(), sv_idx) = 0.0;
    }
  }

  // Set statevector object back to initial state
  statev->update_state(initial_sv);

  return Spectrum(rad.spectral_domain(), 
		  SpectralRange(res, rad.spectral_range().units()));
}
