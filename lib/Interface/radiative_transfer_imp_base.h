#ifndef RADIATIVE_TRANSFER_IMP_BASE_H
#define RADIATIVE_TRANSFER_IMP_BASE_H
#include "radiative_transfer_retrievable.h"
#include "sub_state_vector_array.h"

namespace FullPhysics {
/****************************************************************//**
  As a design principle, we have each base class with the absolutely
  minimum interface needed for use from the rest of the system. This
  allows us to support any future code that supports this minimum 
  interface.
  
  However, almost always you will want to derive from this class 
  instead. See PressureImpBase for a more complete discussion of this.
*******************************************************************/
class RadiativeTransferImpBase: public SubStateVectorArray<RadiativeTransferRetrievable> {
public:
  virtual ~RadiativeTransferImpBase() {}
  virtual boost::shared_ptr<RadiativeTransferRetrievable> clone() const = 0;
 
//-----------------------------------------------------------------------
/// For the sake of being able to return a Spectrum class from Python
/// The reflectance_ptr method here serves the purpose that the radiance
/// method normally would. The reflectance method in this implementation
/// simply calls the reflectance_ptr method for doing the actual work.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Spectrum> reflectance_ptr(const SpectralDomain& Spec_domain, int Spec_index, 
						   bool Skip_jacobian = false) const = 0;

  virtual Spectrum reflectance(const SpectralDomain& Spec_domain, int Spec_index, 
			    bool Skip_jacobian = false) const
  {
    Spectrum res = *reflectance_ptr(Spec_domain, Spec_index, Skip_jacobian);
    return res;
  }

//-----------------------------------------------------------------------
/// Print to stream. The default calls the function "desc" that returns
/// a string. This gives cleaner interface for deriving from this class
/// in python, but most C++ classes will want to override this function
/// rather than using desc.
//-----------------------------------------------------------------------
  virtual void print(std::ostream& Os, bool Short_form = false) const { Os << desc(); }

//-----------------------------------------------------------------------
/// Description of object, to be printed to stream. This gives a cleaner
/// interface for deriving from python.
//-----------------------------------------------------------------------
  virtual std::string desc() const { return "RadiativeTransferImpBase"; }

protected:

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------
  RadiativeTransferImpBase() {}

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() and used_flag() values.
//-----------------------------------------------------------------------
  RadiativeTransferImpBase(const blitz::Array<double, 1>& Coeff, 
			   const blitz::Array<bool, 1>& Used_flag)
    : SubStateVectorArray<RadiativeTransferRetrievable>(Coeff, Used_flag) {}

private:
};
}
#endif
