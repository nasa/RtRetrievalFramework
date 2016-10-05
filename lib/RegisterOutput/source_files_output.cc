#include "source_files_output.h"
#include "fp_exception.h"
// This won't be needed once version 3 becomes the default for boost filesystem
#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(SourceFilesOutput, RegisterOutputBase)
.def(luabind::constructor<const std::string&, const std::vector<std::string>&, const std::vector<std::string>&>())
REGISTER_LUA_END()
#endif

void SourceFilesOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  if(dataset_names_.size() != file_names_.size())
    throw Exception("dataset names and file name vectors need to be the same size");

  for(int name_idx = 0; name_idx < (int) dataset_names_.size(); name_idx++) {
    // Store basename of file, not full path
    std::string base_name = boost::filesystem::path(file_names_[name_idx]).filename().string();
    out->register_data_source(group_name_ + "/" + dataset_names_[name_idx], base_name);
  }
}
