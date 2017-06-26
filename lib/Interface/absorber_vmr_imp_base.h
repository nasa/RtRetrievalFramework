#ifndef ABSORBER_VMR_IMP_BASE_H
#define ABSORBER_VMR_IMP_BASE_H
#include "absorber_vmr.h"
#include "sub_state_vector_array.h"
#include <boost/function.hpp>

namespace FullPhysics {
/****************************************************************//**
  As a design principle, we have each base class with the absolutely
  minimum interface needed for use from the rest of the system. This
  allows us to support any future code that supports this minimum 
  interface.
  
  However, almost always you will want to derive from this class 
  instead. See PressureImpBase for a more complete discussion of this.
*******************************************************************/
class AbsorberVmrImpBase: public SubStateVectorArray<AbsorberVmr> {
public:
  virtual ~AbsorberVmrImpBase() {}
  virtual std::string gas_name() const {return gas_name_;}
  virtual AutoDerivative<double> 
  volume_mixing_ratio(const AutoDerivative<double>& P) const
  { fill_cache(); return vmr(P); }
  virtual boost::shared_ptr<AbsorberVmr> clone() const = 0;
  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const = 0;
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
  virtual std::string desc() const { return "AbsorberVmrImpBase"; }

  virtual blitz::Array<bool, 1> state_used() const
  {
    blitz::Array<bool, 1> res(sv_full.rows());
    res = false;
    if(state_vector_start_index() < res.rows())
      res(blitz::Range(state_vector_start_index(), 
           state_vector_start_index() + sub_vector_size() - 1)) = true;
    return res;
  }
protected:
//-----------------------------------------------------------------------
/// If this is true, the recalculate the vmr the next time we
/// need it.
//-----------------------------------------------------------------------
  mutable bool cache_stale;

//-----------------------------------------------------------------------
/// The cached volumn mixing ration. This should be filled in by
/// derived classes when calc_vmr() is called.
//-----------------------------------------------------------------------
  mutable boost::function<AutoDerivative<double>(AutoDerivative<double>)> vmr;

//-----------------------------------------------------------------------
/// Derived classes should provide a function to fill in vmr when this is 
/// called.
//-----------------------------------------------------------------------
  virtual void calc_vmr() const = 0;

//-----------------------------------------------------------------------
/// Initialize object.
//-----------------------------------------------------------------------

  void init(const std::string Gas_name,
	    const blitz::Array<double, 1>& Coeff, 
	    const blitz::Array<bool, 1>& Used_flag,
	    const boost::shared_ptr<Pressure>& Press,
	    bool Mark_according_to_press = true,
	    int Pdep_start = 0)
  { SubStateVectorArray<AbsorberVmr>::init(Coeff, Used_flag, Press,
					   Mark_according_to_press,
					   Pdep_start);
    gas_name_ = Gas_name;
  }

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------

  AbsorberVmrImpBase() : cache_stale(true) { }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() and used_flag() values.
/// See SubStateVectorArray for a discussion of Mark_according_to_press and
/// Pdep_start.
//-----------------------------------------------------------------------
  AbsorberVmrImpBase(const std::string& Gas_name,
		     const blitz::Array<double, 1>& Coeff, 
		     const blitz::Array<bool, 1>& Used_flag,
		     const boost::shared_ptr<Pressure>& Press,
		     bool Mark_according_to_press = true,
		      int Pdep_start = 0)
    : SubStateVectorArray<AbsorberVmr>(Coeff, Used_flag, Press,
				       Mark_according_to_press, Pdep_start),
      cache_stale(true), gas_name_(Gas_name) { }
private:
  void fill_cache() const
  {
    if(cache_stale)
      calc_vmr();
    cache_stale = false;
  }
  std::string gas_name_;
};
}
#endif
