#ifndef AEROSOL_PARAM_OUTPUT_H
#define AEROSOL_PARAM_OUTPUT_H
#include "register_output_base.h"
#include "aerosol_extinction_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the coefficients used for classes that
  inherits from AerosolExtinctionBaseImp
*******************************************************************/
class AerosolParamOutput : public RegisterOutputBase {
public:
//-----------------------------------------------------------------------
/// Constructor. You can optionally pass in the aerosol name to use in
/// the field names, the default is to use A->aerosol_name().
//-----------------------------------------------------------------------
  AerosolParamOutput(const boost::shared_ptr<AerosolExtinctionImpBase>& A,
		     const std::string& Aname = "") : a(A), aname(Aname) {}
  virtual ~AerosolParamOutput() {}
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<AerosolExtinctionImpBase> a;
  std::string aname;
};
}
#endif
