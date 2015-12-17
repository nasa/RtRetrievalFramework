// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "hdf_file_generating.h"
%}
%base_import(generic_object)
%import "hdf_file.i"

%fp_shared_ptr(FullPhysics::HdfFileGenerating);

namespace FullPhysics {
class HdfFileGenerating : public GenericObject {
public:
  HdfFileGenerating(const std::string& Fname);
  void close();
  void abandon();
  HdfFile& hdf_file();
  std::string print_to_string();
};

}

