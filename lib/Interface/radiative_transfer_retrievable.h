#ifndef RADIATIVE_TRANSFER_RETRIEVABLE_H
#define RADIATIVE_TRANSFER_RETRIEVABLE_H

#include "state_vector.h"
#include "observer.h"
#include "radiative_transfer.h"

namespace FullPhysics {
/****************************************************************//**
 Interface class for radiative transfer implementations that
 happen to have retrievable parameters.								
*******************************************************************/

class RadiativeTransferRetrievable : public RadiativeTransfer,
				     virtual public StateVectorObserver,
				     public Observable<RadiativeTransferRetrievable> {
public:
  virtual ~RadiativeTransferRetrievable() {}

  virtual void add_observer(Observer<RadiativeTransferRetrievable>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<RadiativeTransferRetrievable>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os, bool Short_form = false) const 
  { Os << "RadiativeTransferRetrievable";}

};

}
#endif
