#ifndef SWIG_TO_PYTHON_H
#define SWIG_TO_PYTHON_H
#include "swig_type_mapper_base.h"

// Define this to get diagnostic messages printed out
//#define SWIG_TYPE_MAPPER_DIAGNOSTIC
namespace FullPhysics {
//-----------------------------------------------------------------------
/// Function to map from a shared point to a python object.
//-----------------------------------------------------------------------

template<typename T> inline PyObject* 
swig_to_python(const boost::shared_ptr<T>& V)
{

//-----------------------------------------------------------------------
// If pointer is Null, return None.
//-----------------------------------------------------------------------

  if(!V)
    return Py_None;

//-----------------------------------------------------------------------
// If underlying object is a python object wrapped in a
// Swig::Director, return the underlying python object
//-----------------------------------------------------------------------

  Swig::Director* d = dynamic_cast<Swig::Director*>(V.get());
  if(d) {
#ifdef SWIG_TYPE_MAPPER_DIAGNOSTIC    
    std::cerr << "Return underlying python object\n";
#endif
    return d->swig_get_self();
  }

//-----------------------------------------------------------------------
// See if underlying type is registered in swig_type_map. If so, return the
// underlying type
//-----------------------------------------------------------------------

  T& t(*V.get());
  type_index tid(typeid(t));
#ifdef SWIG_TYPE_MAPPER_DIAGNOSTIC    
  std::cerr << tid.name() << "\n";
#endif
  if(swig_type_map.count(tid) != 0) {
#ifdef SWIG_TYPE_MAPPER_DIAGNOSTIC    
    std::cerr << "Trying to_python for " << tid.name() << "\n";
#endif
    return swig_type_map[tid]->to_python(V);
  }

//-----------------------------------------------------------------------
// Otherwise, fall back to returning the type T.
//-----------------------------------------------------------------------

#ifdef SWIG_TYPE_MAPPER_DIAGNOSTIC    
  std::cerr << "Returning most general type\n";
#endif
  return swig_type_map[typeid(T)]->to_python(V);
}

inline PyObject* 
swig_to_python_or_none(const boost::shared_ptr<GenericObject>& V)
{
//-----------------------------------------------------------------------
// If pointer is Null, return None.
//-----------------------------------------------------------------------

  if(!V)
    return Py_None;

//-----------------------------------------------------------------------
// If underlying object is a python object wrapped in a
// Swig::Director, return the underlying python object
//-----------------------------------------------------------------------

  Swig::Director* d = dynamic_cast<Swig::Director*>(V.get());
  if(d) {
#ifdef SWIG_TYPE_MAPPER_DIAGNOSTIC    
    std::cerr << "Return underlying python object\n";
#endif
    return d->swig_get_self();
  }

//-----------------------------------------------------------------------
// See if underlying type is registered in swig_type_map. If so, return the
// underlying type
//-----------------------------------------------------------------------

  GenericObject& t(*V.get());
  type_index tid(typeid(t));
#ifdef SWIG_TYPE_MAPPER_DIAGNOSTIC    
  std::cerr << tid.name() << "\n";
#endif
  if(swig_type_map.count(tid) != 0) {
#ifdef SWIG_TYPE_MAPPER_DIAGNOSTIC    
    std::cerr << "Trying to_python for " << tid.name() << "\n";
#endif
    return swig_type_map[tid]->to_python(V);
  }

//-----------------------------------------------------------------------
// Otherwise, return Py_None
//-----------------------------------------------------------------------

#ifdef SWIG_TYPE_MAPPER_DIAGNOSTIC    
  std::cerr << "Returning None\n";
#endif
  return Py_None;
}

template<typename T> inline PyObject* 
swig_to_python(const boost::shared_ptr<T>* V)
{ return swig_to_python(*V); }
}
#endif
