#ifndef STOKES_COEFFICIENT_IMP_BASE_H
#define STOKES_COEFFICIENT_IMP_BASE_H
#include "stokes_coefficient.h"
#include "sub_state_vector_array.h"

namespace FullPhysics {
/****************************************************************//**
  As a design principle, we have each base class with the absolutely
  minimum interface needed for use from the rest of the system. This
  allows us to support any future code that supports this minimum 
  interface.
  
  However, almost always you will want to derive from this class 
  instead. This provides additional functionality that you will almost
  always want. We support have a subset set of the full StateVector
  provide coefficients for this class, as well as caching the 
  calculation of the stokes coefficient so we only do the calculation when
  something has changed (e.g., the StateVector).

  The more general StokesCoefficient interface can be used to support unusual
  cases that don't match this implementation, for example wrapping 
  an existing third party library that doesn't mesh with how
  PressureImpBase sets things up. But most of the time, you'll want
  to derive from this class.
*******************************************************************/
class StokesCoefficientImpBase: public SubStateVectorArray<StokesCoefficient> {
public:
  virtual ~StokesCoefficientImpBase() {}
  virtual ArrayAd<double, 2> stokes_coefficient() const
  { fill_cache(); return stokes_coeff; }
  virtual boost::shared_ptr<StokesCoefficient> clone() const = 0;
  virtual void update_sub_state_hook() 
  { cache_stale = true; }

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
  virtual std::string desc() const { return "StokesCoefficientImpBase"; }
protected:
//-----------------------------------------------------------------------
/// If this is true, the recalculate the stokes_coeff the next time we
/// need it.
//-----------------------------------------------------------------------
  mutable bool cache_stale;

//-----------------------------------------------------------------------
/// The cached stokes coefficient. This should be filled in by derived
/// classes when calc_stokes_coeff() is called.
//-----------------------------------------------------------------------
  mutable ArrayAd<double, 2> stokes_coeff;

//-----------------------------------------------------------------------
/// Derived classes should provide a function to fill in pgrid when this is 
/// called.
//-----------------------------------------------------------------------
  virtual void calc_stokes_coeff() const = 0;

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------

  StokesCoefficientImpBase() : cache_stale(true) { }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() and used_flag() values.
//-----------------------------------------------------------------------
  StokesCoefficientImpBase(const blitz::Array<double, 1>& Coeff, 
			   const blitz::Array<bool, 1>& Used_flag)
    : SubStateVectorArray<StokesCoefficient>(Coeff, Used_flag),
      cache_stale(true) {}
private:
  void fill_cache() const
  {
    if(cache_stale) {
      stokes_coeff.resize_number_variable(coeff.number_variable());
      calc_stokes_coeff();
    }
    cache_stale = false;
  }
};

}
#endif
