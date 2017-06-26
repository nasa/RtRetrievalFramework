#include "output.h"
#include <boost/foreach.hpp>

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Handle a single type, this is a helper for pass_to_write.
//-----------------------------------------------------------------------

template<class T> 
void Output::pass_to_write_t(const std::string& Dataset_name, boost::any* D)
{
  typedef typename boost::function<T ()> T2;
  if(boost::any_cast<T2>(D)) {
    T2* t = boost::any_cast<T2>(D);
    write_data(Dataset_name, (*t)());
  }
}

//-----------------------------------------------------------------------
/// Go through the various supported types, and write out the data.
//-----------------------------------------------------------------------

void Output::pass_to_write(const std::string& Dataset_name, boost::any* D)
{
//-----------------------------------------------------------------------
/// Go through each type. Only one of these will actually write and
/// output. 
//-----------------------------------------------------------------------

  pass_to_write_t<int>(Dataset_name, D);
  pass_to_write_t<std::string>(Dataset_name, D);
  pass_to_write_t<const char*>(Dataset_name, D);
  pass_to_write_t<int64_t>(Dataset_name, D);
  pass_to_write_t<double>(Dataset_name, D);

  pass_to_write_t<Array<int, 1> >(Dataset_name, D);
  pass_to_write_t<Array<std::string, 1> >(Dataset_name, D);
  pass_to_write_t<Array<const char*, 1> >(Dataset_name, D);
  pass_to_write_t<Array<double, 1> >(Dataset_name, D);

  pass_to_write_t<Array<int, 2> >(Dataset_name, D);
  pass_to_write_t<Array<std::string, 2> >(Dataset_name, D);
  pass_to_write_t<Array<const char*, 2> >(Dataset_name, D);
  pass_to_write_t<Array<double, 2> >(Dataset_name, D);

  pass_to_write_t<Array<int, 3> >(Dataset_name, D);
  pass_to_write_t<Array<std::string, 3> >(Dataset_name, D);
  pass_to_write_t<Array<const char*, 3> >(Dataset_name, D);
  pass_to_write_t<Array<double, 3> >(Dataset_name, D);
}

//-----------------------------------------------------------------------
/// Write out file. 
///
/// A write is intended to be atomic - either it completely succeeds or
/// the output should be cleaned up and nothing produced. This prevents
/// a file with "missing fields" from being generated.
//-----------------------------------------------------------------------

void Output::write()
{
  try {
    start_write();
    typedef std::map<std::string, boost::any>::value_type vtype;
    BOOST_FOREACH(vtype& p, func)
      pass_to_write(p.first, &p.second);
    end_write();
  } catch(...) {
    try {
      end_because_of_error();	// Notify output to clean itself up.
    } catch(...) {
    }
    throw;
  }
}

//-----------------------------------------------------------------------
/// Write out the file, making a best attempt but ignoring all
/// errors. This may result in a file that is missing fields that had
/// an error when writing or collecting the data. This is intended for
/// use by an error dump, when we want to get whatever we can.
//-----------------------------------------------------------------------

void Output::write_best_attempt()
{
  try {
    start_write();
    typedef std::map<std::string, boost::any>::value_type vtype;
    BOOST_FOREACH(vtype& p, func)
      try {
	pass_to_write(p.first, &p.second);
      } catch(...) {		// Ignore all errors
      }
    end_write();
  } catch(...) {		// Ignore all errors.
  }
}
