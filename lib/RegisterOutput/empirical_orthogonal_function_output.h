#ifndef EMPIRICAL_ORTHOGONAL_FUNCTION_OUTPUT_H
#define EMPIRICAL_ORTHOGONAL_FUNCTION_OUTPUT_H
#include "register_output_base.h"
#include "empirical_orthogonal_function.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the EmpiricalOrthogonalFunction class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the EmpiricalOrthogonalFunction class.
*******************************************************************/
class EmpiricalOrthogonalFunctionOutput : public RegisterOutputBase {
public:
  EmpiricalOrthogonalFunctionOutput
  (const boost::shared_ptr<EmpiricalOrthogonalFunction>& E,
   const std::string& Hdf_band_name) 
    : eof(E), 
      hdf_band_name(Hdf_band_name) {}
  virtual ~EmpiricalOrthogonalFunctionOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<EmpiricalOrthogonalFunction> eof;
  std::string hdf_band_name;
};
}
#endif
