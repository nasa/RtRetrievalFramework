#ifndef ABSORBER_ABSCO_OUTPUT_H
#define ABSORBER_ABSCO_OUTPUT_H
#include "register_output_base.h"
#include "absorber_absco.h"
#include "spectral_bound.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the Absorber class that should be
  written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the Absorber class.
*******************************************************************/
class AbsorberAbscoOutput : public RegisterOutputBase {
public:
  AbsorberAbscoOutput(const boost::shared_ptr<AbsorberAbsco>& A,
		      const SpectralBound& Spec_bound)
    : a(A), sb(Spec_bound), num_level(A->pressure().max_number_level()) {}
  virtual ~AbsorberAbscoOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<AbsorberAbsco> a;
  SpectralBound sb;
  int num_level;
};
}
#endif
