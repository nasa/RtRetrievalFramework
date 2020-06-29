#ifndef ILS_TABLE_LINEAR_OUTPUT_H
#define ILS_TABLE_LINEAR_OUTPUT_H
#include "register_output_base.h"
#include "ils_table.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the IlsTableLinear class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the IlsTableLinear class.
*******************************************************************/
class IlsTableLinearOutput : public RegisterOutputBase {
public:
  IlsTableLinearOutput(const boost::shared_ptr<IlsTableLinear>& D,
                    const std::string& Hdf_band_name) 
    : d(D), 
      hdf_band_name(Hdf_band_name) {}
  virtual ~IlsTableLinearOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<IlsTableLinear> d;
  std::string hdf_band_name;
};
}
#endif
