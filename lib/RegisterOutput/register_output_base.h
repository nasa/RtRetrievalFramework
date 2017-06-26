#ifndef REGISTER_OUTPUT_H
#define REGISTER_OUTPUT_H
#include "printable.h"
#include "output.h"

namespace FullPhysics {
/****************************************************************//**
  As described in the Output class, we have a decentralized model of
  producing output for L2 Full Physics. Rather than directly writing a
  value to a file, we register functions that can supply the values
  when requested. This has the advantage of keeping everything in sync
  - we don't have values written out from earlier iterations of the
  algorithm. 

  There isn't anything particularly special to register a function,
  just a call to Output::register_data_source. However as a matter of
  convention we collect all the functions that register data into
  classes found in the RegisterOutput directory, and derived from this
  RegisterOutputBase class. This class doesn't really give any special
  functionality, rather deriving from this class is a statement of
  intent that derived classes what to register output.

  A reasonable implementation would be to have classes that supply
  output data register there intent themselves, so for example
  Atmosphere could have register_output function. However this would
  then couple these classes with our particular output model. You
  could imagine reusing Atmosphere class in other contexts which do
  not want to use this particular output model. So again as a matter
  of convention we use a separate class to register output, in this
  case AtmosphereOutput.  We may decide at some point that these extra
  classes are unnecessarily complicated design, but for now we'll keep
  this division.

  For many classes, we output different information for apriori state
  vector value vs. the final state vector value. So the registration
  is separated as two functions. Alternatively, we could have just had
  two different registration classes, but this is the way we've
  chosen. 

  Note that by convention that we "freeze" the state of the class 
  when we register the apriori_output. This allows for things like the
  StateVector to be changed after wards without changing the apriori
  state. 
*******************************************************************/
class RegisterOutputBase : public Printable<RegisterOutputBase> {
public:
  virtual ~RegisterOutputBase() {}

//-----------------------------------------------------------------------
/// Register portions of class that will be written to output. This is
/// for the final statevector, for classes where apriori and final are
/// written out different.
//-----------------------------------------------------------------------

  virtual void register_output(const boost::shared_ptr<Output>& out) const = 0;

//-----------------------------------------------------------------------
/// Register apriori portions of class. The default is not to have
/// anything written out, but derived classes can override this.
///
/// Note that by convention that we "freeze" the state of the class 
/// when we register the apriori_output. This allows for things like the
/// StateVector to be changed after wards without changing the apriori
/// state. 
//-----------------------------------------------------------------------

  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const { }

//-----------------------------------------------------------------------
/// Print to stream. The default calls the function "desc" that returns
/// a string. This gives cleaner interface for deriving from this class
/// in python, but most C++ classes will want to override this function
/// rather than using desc.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const { Os << desc();}

//-----------------------------------------------------------------------
/// Description of object, to be printed to stream. This gives a cleaner
/// interface for deriving from python.
//-----------------------------------------------------------------------

  virtual std::string desc() const { return "RegisterOutputBase"; }
};
}
  
#endif
