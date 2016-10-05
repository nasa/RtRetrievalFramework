#ifndef FTS_RUN_LOG_OUTPUT_H
#define FTS_RUN_LOG_OUTPUT_H
#include "register_output_base.h"
#include "fts_run_log.h"

namespace FullPhysics {
/****************************************************************//**
  This writes out data from the FtsRunLog to the output.
*******************************************************************/
class FtsRunLogOutput : public RegisterOutputBase {
public:
  FtsRunLogOutput(const std::vector<FtsRunLogRecord>& Run_log) : 
    run_log(Run_log) {}
  virtual ~FtsRunLogOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
  std::vector<FtsRunLogRecord> run_log;
};
}
#endif
