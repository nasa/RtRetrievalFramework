#include "log_timing.h"
#include "logger.h"
#include "fp_exception.h"
#include <boost/regex.hpp>
#include <sys/resource.h>
#include <boost/lexical_cast.hpp>
#include <unistd.h>
#include <fstream>

using namespace FullPhysics;

//-----------------------------------------------------------------------
/// This is the reference CPU speed that we normalize against. We
/// selected the speed of the scf-srv to use as a reference. 
//-----------------------------------------------------------------------

const double LogTiming::cpu_reference_speed = 2992.497;

// See base class for description
void LogTiming::notify_update(const ConnorSolver& Solver)
{
  number_divergent = Solver.number_divergent();
  number_iteration = Solver.number_iteration();
  write_to_log("Iteration");
}

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

LogTiming::LogTiming() 
  : number_divergent(-1), number_iteration(-1) 
{
  const int host_name_max = 1000;
  char hostnamebuf[host_name_max];
  int status = gethostname(hostnamebuf, host_name_max);
  if(status != 0)
    throw Exception("Call to gethostname failed");
  hostname = std::string(hostnamebuf);
  // The /proc/cpuinfo is Linux specific. If we can't find this file,
  // don't worry about it. We'll just use default values for
  // everything.
  cpu_model = "";
  speed_mhz = -1;
  try {
    std::ifstream proc("/proc/cpuinfo");
    if(proc.good()) {
      proc.exceptions(std::ifstream::badbit | std::ifstream::failbit | 
		      std::ifstream::eofbit);
      std::string ln = "";
      boost::smatch m;
      while(!boost::regex_match(ln, m, 
				boost::regex("model name\\s*:\\s*(.*)")))
	getline(proc, ln);
      cpu_model = m[1];
      while(!boost::regex_match(ln, m, boost::regex("cpu MHz\\s*:\\s*(.*)")))
	getline(proc, ln);
      speed_mhz = boost::lexical_cast<double>(m[1]);
    }
  } catch(std::ifstream::failure e) {
    // Clear any partial initialization.
    cpu_model = "";
    speed_mhz = -1;
  }
  if(speed_mhz < 0)
    cpu_factor = 1.0;		// Default to not doing CPU correction
  else 
    cpu_factor = speed_mhz / cpu_reference_speed;
}

//-----------------------------------------------------------------------
/// Write data to disk.
//-----------------------------------------------------------------------

void LogTiming::write_to_log(const std::string& Prefix) const
{
  // Grab all the timing information at once, so we don't have
  // inconsistent times in the report.
  double wtime = wall_clock.elapsed();
  struct rusage r_usage;
  int status = getrusage(RUSAGE_SELF, &r_usage);
  if(status != 0)
    throw Exception("Call to getrusage failed");
  double utime = r_usage.ru_utime.tv_sec + r_usage.ru_utime.tv_usec * 1e-6;
  double stime = r_usage.ru_stime.tv_sec + r_usage.ru_stime.tv_usec * 1e-6;
  
  std::ostringstream log;
  log << "\n---------------------------------------------------------\n"
      << Prefix << " Timing\n";
  log <<   "Host:                    " << hostname << "\n";
  if(speed_mhz > 0)
    log << "CPU Model:               " << cpu_model << "\n"
	<< "Speed:                   " << speed_mhz << " MHz\n"
	<< "CPU Factor:              " << cpu_factor << "\n";
  else {
    log << "Couldn't get CPU model and speed. Assuming CPU Factor of 1\n"
	<< "CPU Factor:              " << cpu_factor << "\n";
  }
  log <<   "Total Wall clock time:   " << wtime << " seconds\n"
      <<   "Total User time:         " << utime << " seconds\n"
      <<   "Total System time:       " << stime << " seconds\n";
  // /proc/self/status is Linux specific. If we can't find this file,
  // don't worry. We just write a message saying this isn't available.
  try {
    std::ifstream proc("/proc/self/status");
      proc.exceptions(std::ifstream::badbit | std::ifstream::failbit | 
		      std::ifstream::eofbit);
      std::string ln = "";
      boost::smatch m;
      while(!boost::regex_match(ln, m,
				boost::regex("VmPeak:\\s*(\\d+)\\s*kB")))
	getline(proc, ln);
      double t = boost::lexical_cast<double>(m[1]);
      log << "Virtual memory maximum:  " << (t / 1024) << " MB\n";
      while(!boost::regex_match(ln, m,
				boost::regex("VmSize:\\s*(\\d+)\\s*kB")))
	getline(proc, ln);
      t = boost::lexical_cast<double>(m[1]);
      log << "Virtual memory current:  " << (t / 1024) << " MB\n";
      while(!boost::regex_match(ln, m,
				boost::regex("VmHWM:\\s*(\\d+)\\s*kB")))
	getline(proc, ln);
      t = boost::lexical_cast<double>(m[1]);
      log << "Resident memory maximum: " << (t / 1024) << " MB\n";
      while(!boost::regex_match(ln, m,
				boost::regex("VmRSS:\\s*(\\d+)\\s*kB")))
	getline(proc, ln);
      t = boost::lexical_cast<double>(m[1]);
      log << "Resident memory current: " << (t / 1024) << " MB\n";
  } catch(std::ifstream::failure e) {
    log << "Memory usage information isn't available on this platform\n"
	<< "(Available on Linux only)\n";
  } 
  if(number_iteration > 0)
    log << "Number successful steps: " << number_iteration << "\n"
	<< "Number divergent steps:  " << number_divergent << "\n"
	<< "Total iteration:         " << number_iteration + number_divergent
	<< "\n"
	<< "Normalized wall clock time per iteration: " 
	<< wtime / (number_iteration + number_divergent) * cpu_factor
	<< " seconds\n"
	<< "Normalized user time per iteration:       " 
	<< utime / (number_iteration + number_divergent) * cpu_factor
	<< " seconds\n"
	<< "Normalized system time per iteration:     " 
	<< stime / (number_iteration + number_divergent) * cpu_factor
	<< " seconds\n";
  log << "---------------------------------------------------------\n";
  Logger::info() << log.str();
}

