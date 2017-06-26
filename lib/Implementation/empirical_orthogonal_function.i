// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "empirical_orthogonal_function.h"
%}
%base_import(sub_state_vector_array)
%base_import(instrument_correction)
%import "array_with_unit.i"
%import "dispersion.i"
%import "hdf_file.i"
%fp_shared_ptr(FullPhysics::EmpiricalOrthogonalFunction);

namespace FullPhysics {
class EmpiricalOrthogonalFunction: 
  public SubStateVectorArray<InstrumentCorrection> {
public:
  EmpiricalOrthogonalFunction(double Coeff, 
			      bool Used_flag,
			      const ArrayWithUnit<double, 1>& Eof_waveform,
			      int Order,
			      const std::string& Band_name,
			      const std::string& Hdf_group = "N/A");
  EmpiricalOrthogonalFunction(double Coeff, 
			      bool Used_flag,
			      const Dispersion& Disp,
			      const HdfFile& Hdf_static_input,
			      int Spec_index,
			      int Sounding_number,
			      int Order,
			      const std::string& Band_name,
			      const std::string& Hdf_group = 
			      "Instrument/EmpiricalOrthogonalFunction_1");
  EmpiricalOrthogonalFunction(double Coeff, 
			      bool Used_flag,
			      const HdfFile& Hdf_static_input,
			      int Spec_index,
			      int Sounding_number,
			      int Order,
			      const std::string& Band_name,
			      const std::string& Hdf_group = 
			      "Instrument/EmpiricalOrthogonalFunction_1");
  EmpiricalOrthogonalFunction(double Coeff, 
			      bool Used_flag,
			      const HdfFile& Hdf_static_input,
			      const ArrayWithUnit<double, 1>& Uncertainty,
			      int Spec_index,
			      int Sounding_number,
			      int Order,
			      const std::string& Band_name,
			      const std::string& Hdf_group = 
			      "Instrument/EmpiricalOrthogonalFunction",
			      double Scale_to_stddev = 1e19);
  virtual std::string state_vector_name_i(int i) const;
  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const;
  %python_attribute(eof, ArrayWithUnit<double, 1>)
  %python_attribute(order, int)
};
}
