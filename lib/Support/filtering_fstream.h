#ifndef FILTERING_FSTREAM_H
#define FILTERING_FSTREAM_H
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>

namespace FullPhysics {
/****************************************************************//**
   This allows filtering stream to created on top of an 
   underlying ofstream. See boost::iostream documentation for details
   on how to create these. This class is intended to be used as a base
   class for other classes that assemble the proper filters.
*******************************************************************/

class FilteringOstream : public std::ostream {
public:
  virtual ~FilteringOstream() {}
protected:

//-----------------------------------------------------------------------
/// Constructor. This opens the given file and sets up the stream
/// buffer, but nothing is placed in the filter. The derived class should
/// set up the streambuf.
//-----------------------------------------------------------------------

  FilteringOstream(const std::string& Fname) 
    : std::ostream(&sb), f(Fname.c_str()) {}

//-----------------------------------------------------------------------
/// Underlying ofstream.
//-----------------------------------------------------------------------

  std::ofstream f;

//-----------------------------------------------------------------------
/// Filtering streambuf.
//-----------------------------------------------------------------------

  boost::iostreams::filtering_ostreambuf sb;
};

/****************************************************************//**
   This allows filtering stream to created on top of an 
   underlying ifstream. See boost::iostream documentation for details
   on how to create these. This class is intended to be used as a base
   class for other classes that assemble the proper filters.
*******************************************************************/

class FilteringIstream : public std::istream {
public:
  virtual ~FilteringIstream() {}
protected:

//-----------------------------------------------------------------------
/// Constructor. This opens the given file and sets up the stream
/// buffer, but nothing is placed in the filter. The derived class should
/// set up the streambuf.
//-----------------------------------------------------------------------

  FilteringIstream(const std::string& Fname) 
    : std::istream(&sb), f(Fname.c_str()) {}

//-----------------------------------------------------------------------
/// Underlying ifstream.
//-----------------------------------------------------------------------

  std::ifstream f;

//-----------------------------------------------------------------------
/// Filtering streambuf.
//-----------------------------------------------------------------------

  boost::iostreams::filtering_istreambuf sb;
};
}
#endif
