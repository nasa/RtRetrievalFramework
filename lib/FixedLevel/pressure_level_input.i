// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "pressure_level_input.h"
%}

%import "heritage_file.i"
%import "hdf_file.i"
%fp_shared_ptr(FullPhysics::PressureLevelInput);

namespace FullPhysics {
class PressureLevelInput {
public:
  PressureLevelInput(const blitz::Array<double, 1>& Press_level);
  PressureLevelInput(const HeritageFile& Heritage_file);
  PressureLevelInput(const HdfFile& Hdf_file, 
		     const std::string& Hdf_group = "Pressure");
  %python_attribute(pressure_level, blitz::Array<double, 1>)
  std::string print_to_string() const;
};
}
