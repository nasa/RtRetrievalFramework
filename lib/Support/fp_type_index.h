#ifndef FP_TYPE_INDEX_H
#define FP_TYPE_INDEX_H
#include <cstring>
#include <typeinfo>
#include <boost/operators.hpp>

namespace FullPhysics {
  /// Dummy class to use as a nice default for type_index
class null_type {
};

/****************************************************************//**
  This is a wrapper around std::type_info that allows it to be used
  as an index in a associative container.

  This is actually in cxx11, but we don't want to depend on using a
  cxx11 compiler. When these become more the standard compiler, we
  can replace this class with std::type_info.
*******************************************************************/
class type_index : public boost::totally_ordered<type_index> {
public:

//-----------------------------------------------------------------------
/// Default constructor.
//-----------------------------------------------------------------------

  type_index() : id(&typeid(null_type)) {}

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

  type_index(std::type_info const& id) : id(&id) {}

//-----------------------------------------------------------------------
/// Comparison operator
//-----------------------------------------------------------------------

  bool operator==(type_index const& other) const 
  { return std::strcmp(id->name(), other.id->name()) == 0;}

//-----------------------------------------------------------------------
/// Comparison operator
//-----------------------------------------------------------------------

  bool operator<(type_index const& other) const 
  {return std::strcmp(id->name(), other.id->name()) < 0;}

//-----------------------------------------------------------------------
/// Return type name.
//-----------------------------------------------------------------------

  std::string name() const {return id->name();}
private:
  std::type_info const* id;
};
}
#endif
