#ifndef LOG_TIMING_H
#define LOG_TIMING_H
#include "printable.h"
#include "connor_solver.h"
#include <boost/timer.hpp>

namespace FullPhysics {

/****************************************************************//**
  This is a helper class that logs basic timing information to the
  Logger.
*******************************************************************/

class LogTiming : public Printable<LogTiming>, public Observer<ConnorSolver> {
public:
  LogTiming();
  virtual ~LogTiming() {};
  virtual void notify_update(const ConnorSolver& Solver);
  virtual void print(std::ostream& Os) { Os << "LogTiming";}
  virtual void write_to_log(const std::string& Prefix = "") const;
private:
  int number_divergent;
  int number_iteration;
  boost::timer wall_clock;
  double cpu_factor;
  double speed_mhz;
  std::string cpu_model;
  std::string hostname;
  static const double cpu_reference_speed;
};
}

#endif
