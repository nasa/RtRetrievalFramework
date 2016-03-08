#include "hdf_sounding_id.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(HdfSoundingId)
.def("frame_number", &HdfSoundingId::frame_number)
.def("sounding_number", &HdfSoundingId::sounding_number)
.def("sounding_id", &HdfSoundingId::sounding_id)
REGISTER_LUA_END()

typedef std::vector<boost::shared_ptr<HdfSoundingId> >::reference 
(std::vector<boost::shared_ptr<HdfSoundingId> >::*vt1)(
        std::vector<boost::shared_ptr<HdfSoundingId> >::size_type);
// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<HdfSoundingId> >::*pbt1)(
        const std::vector<boost::shared_ptr<HdfSoundingId> >::value_type&);
REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<HdfSoundingId> >, VectorHdfSoundingId)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<HdfSoundingId> >::push_back))
.def("size", &std::vector<boost::shared_ptr<HdfSoundingId> >::size)
.def("value", ((vt1) &std::vector<boost::shared_ptr<HdfSoundingId> >::operator[]))
REGISTER_LUA_END()
#endif
