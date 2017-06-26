#ifndef TEMPERATURE_IMP_BASE_H
#define TEMPERATURE_IMP_BASE_H
#include "temperature.h"
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
class TemperatureImpBase: public SubStateVectorArray<Temperature> {
public:
  virtual ~TemperatureImpBase() {}
  virtual AutoDerivativeWithUnit<double> 
  temperature(const AutoDerivativeWithUnit<double>& Press) const
  { fill_cache(); 
    return AutoDerivativeWithUnit<double>(tgrid(Press.convert(units::Pa).value),
					  units::K); 
  }
  virtual boost::shared_ptr<Temperature> clone() const = 0;
  virtual boost::shared_ptr<Temperature> 
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
  virtual std::string desc() const { return "TemperatureImpBase"; }
protected:
//-----------------------------------------------------------------------
/// If this is true, the recalculate the temperature_grid the next time we
/// need it.
//-----------------------------------------------------------------------
  mutable bool cache_stale;

//-----------------------------------------------------------------------
/// The cached temperature grid. This should be filled in by derived classes
/// when calc_temperature_grid() is called. This should map pressure
/// is Pascal to Temperature in Kelvin.
//-----------------------------------------------------------------------
  mutable boost::function<AutoDerivative<double>(AutoDerivative<double>)> tgrid;

//-----------------------------------------------------------------------
/// Derived classes should provide a function to fill in tgrid when this is 
/// called.
//-----------------------------------------------------------------------
  virtual void calc_temperature_grid() const = 0;

//-----------------------------------------------------------------------
/// Initialize object.
//-----------------------------------------------------------------------

  void init(const blitz::Array<double, 1>& Coeff, 
	    const blitz::Array<bool, 1>& Used_flag,
	    const boost::shared_ptr<Pressure>& Press,
	    bool Mark_according_to_press = true,
	    int Pdep_start = 0)
  { SubStateVectorArray<Temperature>::init(Coeff, Used_flag, Press,
					   Mark_according_to_press,
					   Pdep_start);
  }

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------

  TemperatureImpBase() : cache_stale(true) { }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() and used_flag() values.
/// See SubStateVectorArray for a discussion of Mark_according_to_press and
/// Pdep_start.
//-----------------------------------------------------------------------
  TemperatureImpBase(const blitz::Array<double, 1>& Coeff, 
		     const blitz::Array<bool, 1>& Used_flag,
		     const boost::shared_ptr<Pressure>& Press,
		     bool Mark_according_to_press = true,
		     int Pdep_start = 0)
    : SubStateVectorArray<Temperature>(Coeff, Used_flag, Press,
				       Mark_according_to_press, Pdep_start),
      cache_stale(true) { }
private:
  void fill_cache() const
  {
    if(cache_stale)
      calc_temperature_grid();
    cache_stale = false;
  }

};
}
#endif
