#ifndef AEROSOL_PROPERTY_IMP_BASE_H
#define AEROSOL_PROPERTY_IMP_BASE_H
#include "aerosol_property.h"
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
class AerosolPropertyImpBase: public SubStateVectorArray<AerosolProperty> {
public:
  virtual ~AerosolPropertyImpBase() {}
  virtual boost::shared_ptr<AerosolProperty> clone() const = 0;
  virtual boost::shared_ptr<AerosolProperty> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<RelativeHumidity>& Rh) const = 0;
  virtual ArrayAd<double,1> extinction_coefficient_each_layer(double wn) 
    const = 0;
  virtual ArrayAd<double, 1> scattering_coefficient_each_layer(double wn)
    const = 0;
  virtual ArrayAd<double, 3> 
  phase_function_moment_each_layer(double wn, int nmom = -1, 
				   int nscatt = -1) const = 0;

//-----------------------------------------------------------------------
/// Print to stream. The default calls the function "desc" that returns
/// a string. This gives cleaner interface for deriving from this class
/// in python, but most C++ classes will want to override this function
/// rather than using desc.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const { Os << desc(); }

//-----------------------------------------------------------------------
/// Description of object, to be printed to stream. This gives a cleaner
/// interface for deriving from python.
//-----------------------------------------------------------------------

  virtual std::string desc() const { return "AerosolPropertyImpBase"; }

//-----------------------------------------------------------------------
/// Returns the value of the coefficients used to generate the aerosol
/// property
//-----------------------------------------------------------------------

  blitz::Array<double, 1> aerosol_parameter() const
  {
    return coefficient().value();
  }

//-----------------------------------------------------------------------
/// Returns the uncertainty of the aerosol type coefficients
//-----------------------------------------------------------------------

  blitz::Array<double, 1> aerosol_parameter_uncertainty() const
  {
    blitz::Array<double, 1> uncert(coefficient().rows());
    for(int i = 0; i < sv_cov_sub.rows(); i++)
      uncert(i) = (sv_cov_sub(i,i) > 0 ? sqrt(sv_cov_sub(i, i)) : 0.0);
    return uncert;
  }

protected:
  // Don't think we need a cache, but if we end up needing it can add
  // this like we have in AerosolExtinctionImpBase.
  // mutable bool cache_stale;

//-----------------------------------------------------------------------
/// Initialize object.
//-----------------------------------------------------------------------

  void init(const blitz::Array<double, 1>& Coeff, 
	    const blitz::Array<bool, 1>& Used_flag)
  { 
    SubStateVectorArray<AerosolProperty>::init(Coeff, Used_flag);
  }

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------

  AerosolPropertyImpBase() { }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() and used_flag() values.
/// See SubStateVectorArray for a discussion of Mark_according_to_press and
/// Pdep_start.
//-----------------------------------------------------------------------
  AerosolPropertyImpBase(const blitz::Array<double, 1>& Coeff, 
			 const blitz::Array<bool, 1>& Used_flag)
    : SubStateVectorArray<AerosolProperty>(Coeff, Used_flag)
      { }
};
}
#endif
