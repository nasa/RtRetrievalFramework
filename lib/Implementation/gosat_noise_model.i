// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "gosat_noise_model.h"
#include "acos_sounding_id.h"
%}
%base_import(noise_model)
%import "heritage_file.i"
%import "hdf_file.i"
%import "acos_sounding_id.i"
%import "instrument.i"
%fp_shared_ptr(FullPhysics::GosatNoiseModel);

namespace FullPhysics {
class GosatNoiseModel: public NoiseModel {
public:
  GosatNoiseModel(const HeritageFile& Noise_file);
  GosatNoiseModel(const HeritageFile& Noise_file, const HeritageFile& Emp_Coeff_File);
  GosatNoiseModel(const HdfFile& Noise_file, const AcosSoundingId& Sounding_Id, const Instrument& Inst);
  GosatNoiseModel(const HdfFile& Noise_file, const AcosSoundingId& Sounding_Id, const Instrument& Inst, const HeritageFile& Emp_Coeff_File);
  virtual blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const;
};
}
