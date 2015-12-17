#ifndef DISPERSION_FIT_OUTPUT_H
#define DISPERSION_FIT_OUTPUT_H
#include "register_output_base.h"
#include "dispersion_fit.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the DispersionFit class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the DispersionFit class.
*******************************************************************/
class DispersionFitOutput : public RegisterOutputBase {
public:
  DispersionFitOutput(const boost::shared_ptr<DispersionFit>& D) 
    : d(D) {}
  virtual ~DispersionFitOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<DispersionFit> d;
};
}
#endif
