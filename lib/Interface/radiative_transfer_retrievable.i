// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "radiative_transfer_retrievable.h"
#include "sub_state_vector_array.h"
%}

%base_import(radiative_transfer)
%base_import(state_vector)
%base_import(sub_state_vector_array)
%base_import(observer)
%import "spectral_domain.i"

%fp_shared_ptr(FullPhysics::RadiativeTransferRetrievable)

namespace FullPhysics {
  class RadiativeTransferRetrievable;
}

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::RadiativeTransferRetrievable>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::RadiativeTransferRetrievable>)
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::RadiativeTransferRetrievable>);


namespace FullPhysics {
%template(ObservableRadiativeTransferRetrievable) FullPhysics::Observable<RadiativeTransferRetrievable>;
%template(ObserverRadiativeTransferRetrievable) FullPhysics::Observer<RadiativeTransferRetrievable>;

class RadiativeTransferRetrievable : public RadiativeTransfer,
				     virtual public StateVectorObserver,
				     public Observable<RadiativeTransferRetrievable> {
public:
  virtual ~RadiativeTransferRetrievable();

  virtual void add_observer(Observer<RadiativeTransferRetrievable>& Obs);
  virtual void remove_observer(Observer<RadiativeTransferRetrievable>& Obs);
};

%template(SubStateVectorArrayRadiativeTransfer) 
     FullPhysics::SubStateVectorArray<RadiativeTransferRetrievable>;

}
