#ifndef SWIG_TYPE_MAPPER_H
#define SWIG_TYPE_MAPPER_H
#include "swig_type_mapper_base.h"
#include "swig_to_python.h"

namespace FullPhysics {

/****************************************************************//**
  This is the implementation of SwigTypeMapperBase
*******************************************************************/

template<class T> class SwigTypeMapper : public SwigTypeMapperBase {
public:
  SwigTypeMapper(const char* Typename) {
    sinfo = SWIG_TypeQuery(Typename);
#ifdef SWIG_TYPE_MAPPER_DIAGNOSTIC    
    std::cerr << "sinfo - " << sinfo->str << "\n";
#endif
  }
  virtual PyObject* to_python(const boost::shared_ptr<GenericObject>& V) 
  {
    boost::shared_ptr<T> v2 = boost::dynamic_pointer_cast<T>(V);
    boost::shared_ptr<T> *v3 = v2 ? new boost::shared_ptr<T>(v2) : 0;
    return SWIG_NewPointerObj(SWIG_as_voidptr(v3), sinfo, 
			      SWIG_POINTER_OWN);
  }
private:
  swig_type_info* sinfo;
};
}
#endif
