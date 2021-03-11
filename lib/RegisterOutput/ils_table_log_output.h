#ifndef ILS_TABLE_LOG_OUTPUT_H
#define ILS_TABLE_LOG_OUTPUT_H
#include "register_output_base.h"
#include "ils_table.h"

namespace FullPhysics {
/****************************************************************//**
  This registers the portions of the IlsTableLog class that
  should be written as output.

  See the discussion in RegisterOutputBase why this isn't just part of
  the IlsTableLog class.
*******************************************************************/
class IlsTableLogOutput : public RegisterOutputBase {
public:
  IlsTableLogOutput(const boost::shared_ptr<IlsTableLog>& D,
                    const std::string& Hdf_band_name) 
    : d(D), 
      hdf_band_name(Hdf_band_name) {}
  virtual ~IlsTableLogOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
private:
  boost::shared_ptr<IlsTableLog> d;
  std::string hdf_band_name;
};
}
#endif
