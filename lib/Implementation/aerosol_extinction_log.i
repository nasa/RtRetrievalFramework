// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "aerosol_extinction_log.h"
  %}
%base_import(aerosol_extinction_imp_base)
%import "pressure.i"
%fp_shared_ptr(FullPhysics::AerosolExtinctionLog)

// Force to be not abstract, SWIG had troubles seeing that the clone methods ARE implemented below
%feature("notabstract") AerosolExtinctionLog;

namespace FullPhysics {
class AerosolExtinctionLog : public AerosolExtinctionImpBase { 
public:
  AerosolExtinctionLog(const boost::shared_ptr<Pressure>& Press,
			  const blitz::Array<bool, 1>& Flag, 
			  const blitz::Array<double, 1>& Aext,
			  const std::string& Aerosol_name);
  virtual boost::shared_ptr<AerosolExtinction> clone() const;
  virtual boost::shared_ptr<AerosolExtinction> clone
  (const boost::shared_ptr<Pressure>& P) const;
  virtual std::string state_vector_name_i(int i) const;
protected:
  virtual void calc_aerosol_extinction() const;
};
}

