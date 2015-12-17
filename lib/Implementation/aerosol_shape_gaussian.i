// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "aerosol_shape_gaussian.h"
%}
%base_import(aerosol_extinction_imp_base)
%import "pressure.i"
%fp_shared_ptr(FullPhysics::AerosolShapeGaussian)
namespace FullPhysics {
class AerosolShapeGaussian : public AerosolExtinctionImpBase {
public:
  AerosolShapeGaussian(const boost::shared_ptr<Pressure>& Press,
		       const blitz::Array<bool, 1>& Flag, 
		       const blitz::Array<double, 1>& Coeffs,
		       const std::string& Aerosol_name,
		       const bool Linear_AOD);
  virtual boost::shared_ptr<AerosolExtinction> clone() const;
  virtual boost::shared_ptr<AerosolExtinction> clone
  (const boost::shared_ptr<Pressure>& P) const;
  virtual std::string state_vector_name_i(int i) const;
protected:
  virtual void calc_aerosol_extinction() const;
};
}
