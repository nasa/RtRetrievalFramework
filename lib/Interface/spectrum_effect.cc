#include "spectrum_effect.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(SpectrumEffect)
REGISTER_LUA_END()

typedef std::vector<boost::shared_ptr<SpectrumEffect> >::reference 
(std::vector<boost::shared_ptr<SpectrumEffect> >::*vsevt)(std::vector<boost::shared_ptr<SpectrumEffect> >::size_type);

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<SpectrumEffect> >::*pbt1)(
        const std::vector<boost::shared_ptr<SpectrumEffect> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<SpectrumEffect> >, VectorSpectrumEffect)
.def(luabind::constructor<>())
.def("size", &std::vector<boost::shared_ptr<SpectrumEffect> >::size)
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<SpectrumEffect> >::push_back))
.def("value", ((vsevt) &std::vector<boost::shared_ptr<SpectrumEffect> >::operator[]))
REGISTER_LUA_END()

typedef std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::reference 
(std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::*vvsevt)(std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::size_type);

typedef void(std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::*pbt2)(
        const std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >, VectorVectorSpectrumEffect)
.def(luabind::constructor<>())
.def("size", &std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::size)
.def("push_back", ((pbt2) &std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::push_back))
.def("value", ((vvsevt) &std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::operator[]))
REGISTER_LUA_END()

#endif
