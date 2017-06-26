#ifndef OUTPUT_H
#define OUTPUT_H
#include "printable.h"
#include <boost/function.hpp>
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/any.hpp>
#include <map>
#include <string>
#include <stdint.h>

namespace FullPhysics {
// Helper class to handle level padding.
// Don't have Doxygen document this class.
/// @cond
template<class T> class LevelPad {
public:
  LevelPad(const boost::function<blitz::Array<T, 1> ()>& F, T Fill, int Fullsize)
    : f(F), fill(Fill), fullsize(Fullsize) {}
  blitz::Array<T, 1> operator()() const
  {
    blitz::Array<T, 1> res(fullsize);
    res = fill;
    blitz::Array<T, 1> t = f();
    res(blitz::Range(0, t.rows() - 1)) = t;
    return res;
  }
private:
  boost::function<blitz::Array<T, 1> ()> f;
  T fill;
  int fullsize;
};

template<class T> class LevelPad2 {
public:
  LevelPad2(const boost::function<blitz::Array<T, 2> ()>& F, T Fill, int Fullsize)
    : f(F), fill(Fill), fullsize(Fullsize) {}
  blitz::Array<T, 2> operator()() const
  {
    blitz::Array<T, 2> res(fullsize);
    res = fill;
    blitz::Array<T, 2> t = f();
    res(blitz::Range(0, t.rows() - 1),
	blitz::Range(0, t.rows() - 1)) = t;
    return res;
  }
private:
  boost::function<blitz::Array<T, 2> ()> f;
  T fill;
  int fullsize;
};
/// @endcond



/****************************************************************//**
  This is the base class for classes that write output for Level 2
  Full Physics. Specific derived classes are used to write out files
  (e.g., OutputHdf, OutputHeritage), this class just captures the
  common behavior.

  The output is designed to be decentralized. Classes register
  themselves as being able to supply data to be written out when
  requested, and the OutputHeritage class will request this data when
  it wants to write out the state of the system. This "write on
  demand" design allows for scenarios such as "dump data out when an
  error occurs" and "dump data on each iteration of the retrieval
  solver", and "checkpoint results".

  In many cases, the dataset and metadata is supplied either as a
  constant unchanging value, or as a pointer to a member function of
  an object supplied as a boost::shared_ptr<T>. The functions 
  "register_data_source" have been overloaded to handle these common
  cases, most of the time these are the only functions you need. If 
  you want to specify this with a more general boost::function, you 
  can do so also.

  The derived class need to supply various write_xxx functions to act
  as a data sink for this data.

  This class only supports a fixed number of data types (int, double,
  std::string) and array rank (scalar, 1d, 2d, 3d). We also have a
  int64_t scalar, needed for the sounding id. We can easily extend
  this to other types and ranks if needed. For developers, you add a
  new types by adding a new write_data function to this class, along
  with adding the type to function "pass_to_write".

  A write is intended to be atomic - either it completely succeeds or
  the output should be cleaned up and nothing produced. This prevents
  a file with "missing fields" from being generated.

  In the case of a processing error, we may want to make an attempt to
  dump as much of the data as possible to help with diagnostics. A
  separate "write_best_attempt" routines is supplied, this will write
  out whatever we can, ignoring all errors. This can results in
  partial files, but in the case of a diagnostic file whatever we can
  get is better than nothing.
*******************************************************************/

class Output : public Printable<Output> {
public:
  virtual ~Output() {};
  virtual void print(std::ostream& Os) const { Os << "Output";}

//-----------------------------------------------------------------------
/// A common way to supply the metadata source is with a shared_ptr to
/// an object, and a pointer to a member function. This supplies a
/// simplified interface to this. 
//-----------------------------------------------------------------------

  template<class S, class T, int D> void
  register_data_source(const std::string& Dataset_name,
		       blitz::Array<T, D> (S::*Pmf)() const,
		       const boost::shared_ptr<S>& Src)
  {
    boost::function<blitz::Array<T, D> ()> f = boost::bind(Pmf, Src);
    register_data_source(Dataset_name, f);
  }

//-----------------------------------------------------------------------
/// There are several fields that are generated on the active levels
/// only. For output, we want to pad this to the full number of
/// levels. This utility routine will pad the given size using the
/// given fill value.
//-----------------------------------------------------------------------

  template<class S, class T> void
  register_data_source_pad(const std::string& Dataset_name,
			   blitz::Array<T, 1> (S::*Pmf)() const,
			   const boost::shared_ptr<S>& Src,
			   int Full_size,
			   T Fill_value)
  {
    boost::function<blitz::Array<T, 1> ()> f = boost::bind(Pmf, Src);
    boost::function<blitz::Array<T, 1> ()> f2 =
      LevelPad<T>(f, Fill_value, Full_size);
    register_data_source(Dataset_name, f2);
  }
  
  template<class S, class T> void
  register_data_source_pad(const std::string& Dataset_name,
			   blitz::Array<T, 2> (S::*Pmf)() const,
			   const boost::shared_ptr<S>& Src,
			   int Full_size,
			   T Fill_value)
  {
    boost::function<blitz::Array<T, 2> ()> f = boost::bind(Pmf, Src);
    boost::function<blitz::Array<T, 2> ()> f2 =
      LevelPad2<T>(f, Fill_value, Full_size);
    register_data_source(Dataset_name, f2);
  }
  
//-----------------------------------------------------------------------
/// A common way to supply the metadata source is with a shared_ptr to
/// an object, and a pointer to a member function. This supplies a
/// simplified interface to this. 
//-----------------------------------------------------------------------

  template<class S, class T> void
  register_data_source(const std::string& Dataset_name,
		       T (S::*Pmf)() const,
		       const boost::shared_ptr<S>& Src)
  {
    boost::function<T ()> f = boost::bind(Pmf, Src);
    register_data_source(Dataset_name, f);
  }
  
//-----------------------------------------------------------------------
/// Handling for when data is constant.
//-----------------------------------------------------------------------

  template<class T, int D> void
  register_data_source(const std::string& Dataset_name,
		       const blitz::Array<T, D>& Val)
  {
    boost::function<blitz::Array<T, D> ()> f = 
      boost::lambda::constant(Val.copy());
    register_data_source(Dataset_name, f);
  }

//-----------------------------------------------------------------------
/// Handling for when data is constant.
//-----------------------------------------------------------------------

  template<class T> void
  register_data_source(const std::string& Dataset_name,
		       const T& Val)
  {
    boost::function<T ()> f = boost::lambda::constant(Val);
    register_data_source(Dataset_name, f);
  }

//-----------------------------------------------------------------------
/// Most general way to enter a Data source.
//-----------------------------------------------------------------------

  template<class T> 
  void register_data_source(const std::string& Dataset_name,
			    boost::function<T> f)
  { func[Dataset_name] = f;}
  void write();
  void write_best_attempt();
protected:

//-----------------------------------------------------------------------
/// Notify derived class that we are starting to write data. Default
/// is to do nothing, but derived classes can override this.
//-----------------------------------------------------------------------

  virtual void start_write() {};

//-----------------------------------------------------------------------
/// Notify when we are done. Default is to do nothing, but derived
/// classes can override this.
//-----------------------------------------------------------------------

  virtual void end_write() {};

//-----------------------------------------------------------------------
/// Notify when an error occurred. Derived classes should clean up any
/// partially generated output.
//-----------------------------------------------------------------------

  virtual void end_because_of_error() {}

//-----------------------------------------------------------------------
/// Various write functions that derived classes 
//-----------------------------------------------------------------------
  virtual void write_data(const std::string& Dataset_name, int Val) = 0;
  virtual void write_data(const std::string& Dataset_name, int64_t Val) = 0;
  virtual void write_data(const std::string& Dataset_name, double Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
			  const std::string& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
			  const char* Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<int, 1>& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<std::string, 1>& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<const char*, 1>& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<double, 1>& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<int, 2>& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<std::string, 2>& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<const char*, 2>& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<double, 2>& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<int, 3>& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<std::string, 3>& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<const char*, 3>& Val) = 0;
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<double, 3>& Val) = 0;
private:
  template<class T> 
  void pass_to_write_t(const std::string& Dataset_name, boost::any* D);
  void pass_to_write(const std::string& Dataset_name, boost::any* D);
  std::map<std::string, boost::any> func;
};

/****************************************************************//**
  Most of the time the write_data needed by Output is best done
  through a template. However you can't actually have virtual
  templates. This base class implements each of the write_data member
  functions by a call to a function "write_data_t" that should be
  derived in a class T.
*******************************************************************/
template<class T> class OutputTemplate : public Output 
{
public:
  virtual ~OutputTemplate() {}
protected:
  virtual void write_data(const std::string& Dataset_name, int Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, int64_t Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, double Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
			  const std::string& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
			  const char* Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<int, 1>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<std::string, 1>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<const char*, 1>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<double, 1>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<int, 2>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<std::string, 2>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<const char*, 2>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<double, 2>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<int, 3>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<std::string, 3>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<const char*, 3>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
  virtual void write_data(const std::string& Dataset_name, 
		  const blitz::Array<double, 3>& Val)
  { ((T*)this)->write_data_t(Dataset_name, Val); }
};

class OutputDouble : public GenericObject {
public:
  OutputDouble() {}
  virtual ~OutputDouble() {}
  virtual double f() const = 0;
};

class OutputBlitz1d : public GenericObject {
public:
  OutputBlitz1d() {}
  virtual ~OutputBlitz1d() {}
  virtual blitz::Array<double,1> f() const = 0;
};

class OutputBlitz2d  : public GenericObject {
public:
  OutputBlitz2d() {}
  virtual ~OutputBlitz2d() {}
  virtual blitz::Array<double,2> f() const = 0;
};
}
  
#endif
