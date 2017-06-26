#ifndef GAS_VMR_APRIORI_OUTPUT_H
#define GAS_VMR_APRIORI_OUTPUT_H

#include "register_output_base.h"
#include "gas_vmr_apriori.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the GasVmrApriori class that
  should be written as output.
*******************************************************************/
class GasVmrAprioriOutput : public RegisterOutputBase {
public:
  GasVmrAprioriOutput(const boost::shared_ptr<GasVmrApriori>& Gas_apriori) 
    : gas_apriori(Gas_apriori) {}
  virtual ~GasVmrAprioriOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<GasVmrApriori> gas_apriori;
};
}
#endif
