#include "spectrum_effect.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(SpectrumEffect)
REGISTER_LUA_END()

typedef std::vector<boost::shared_ptr<SpectrumEffect> >::reference 
(std::vector<boost::shared_ptr<SpectrumEffect> >::*vsevt)(std::vector<boost::shared_ptr<SpectrumEffect> >::size_type);
REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<SpectrumEffect> >,
			VectorSpectrumEffect)
.def(luabind::constructor<>())
.def("size", &std::vector<boost::shared_ptr<SpectrumEffect> >::size)
.def("push_back", &std::vector<boost::shared_ptr<SpectrumEffect> >::push_back)
.def("value", ((vsevt) &std::vector<boost::shared_ptr<SpectrumEffect> >::operator[]))
REGISTER_LUA_END()

typedef std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::reference 
(std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::*vvsevt)(std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::size_type);
REGISTER_LUA_CLASS_NAME(std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >,
			VectorVectorSpectrumEffect)
.def(luabind::constructor<>())
.def("size", &std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::size)
.def("push_back", &std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::push_back)
.def("value", ((vvsevt) &std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >::operator[]))
REGISTER_LUA_END()

#endif
