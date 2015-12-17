// This is the Level 2 Full Physics main program.
#include "fp_exception.h"
#include "l2_fp_configuration_lua.h"
#include "fp_logger.h"
#include "log_timing.h"
#include <signal.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <fenv.h>
#include <H5Cpp.h>
#include <ctime>
#include "version.h"

using namespace FullPhysics;

static boost::shared_ptr<Output> output_error;
static boost::shared_ptr<LogTiming> log_timing_ptr;
// This is a handler in case we get a floating point exception.
void floating_point_exception(int signum, siginfo_t *sip, void *scp)
{
  std::cerr << "Exception thrown by Full Physics code:\n"
            << "   Floating point exception\n"
	    << "   Address exception occurred: " << sip->si_addr << "\n"
	    << "   Exception type:\n" << "      ";
  int fe_code = sip->si_code;
  switch(fe_code) {
  case FPE_INTDIV:
    std::cerr << "Integer divide by zero\n";
    break;
  case FPE_INTOVF:
    std::cerr << "Integer overflow\n";
    break;
  case FPE_FLTDIV:
    std::cerr << "Floating point divide by zero\n";
    break;
  case FPE_FLTOVF:
    std::cerr << "Floating point overflow\n";
    break;
  case FPE_FLTUND:
    std::cerr << "Floating point underflow\n";
    break;
  case FPE_FLTRES:
    std::cerr << "Floating point inexact result\n";
    break;
  case FPE_FLTINV:
    std::cerr << "Floating point invalid operation\n";
    break;
  case FPE_FLTSUB:
    std::cerr << "Subscript out of range\n";
    break;
  default:
    std::cerr << "Unknown floating point exception.\n";
  }
  if(output_error)
    output_error->write_best_attempt();
  if(log_timing_ptr)
    log_timing_ptr->write_to_log("Final Error");
  exit(5);
}

int main(int Argc, char** Argv)
{
  // Need to create a logger before the one created by the configuration file
  // since there may be items needing to be logged as a result of the actions
  // of configuration heritage
  Logger::set_implementation(boost::shared_ptr<LogImp>
                             (new FpLogger(LogImp::INFO)));
  log_timing_ptr.reset(new LogTiming);
  LogTiming& log_timing = *log_timing_ptr;
  time_t tm;                        // Scratch variable use to get time information
  try {
    bool help = false;
    bool save_test = false;
    bool show_version = false;
    std::string save_test_file;
    char ch;
    while((ch = getopt(Argc, Argv, "t:vh")) != -1) {
      switch(ch) {
      case 't':
        save_test = true;
        save_test_file = optarg;
        break;
      case 'v':
        show_version = true;
        break;
      default:
        help = true;
      }
    }

    Argc -= optind;
    Argv += optind;

    if(help) {
      std::cerr <<
"Usage: l2_fp [-h] [-v] [-t <test file name>] [<configuration file>] \n"
"   [<output file>]\n"
"\n"
"This program is used to do a Level 2 Full Physics retrieval. You\n"
"can optionally specify a configuration file to use and an output\n"
"file. The default is './oco_l2.run' for the configuration file, and\n"
"'./out.h5' for the output file.\n"
"\n"
"The configuration can be the old heritage format, or the new Lua\n"
"config file. If the file name ends with \".lua\" then we use the\n"
"new Lua format, otherwise it is assumed to tbe the old format.\n"
"\n"
"As a special backdoor, you can optionally specify -t. This causes\n"
"the final state of the ConnorSolver object to be saved as the \n"
"given file name. This file can then be used in unit tests. You\n"
"wouldn't normally specify this option.\n";
      return 1;
    }

    // Check for a special Lua file named version.lua somewhere
    // in the LUA_PATH. Use this to extract a version for the 
    // Lua config.
    boost::shared_ptr<LuaState> ls(new LuaState("."));
    std::string require_cmd = "require(\"version\")";
    int status = luaL_dostring(ls->lua_state(), require_cmd.c_str());
    std::string lua_cm_ver = "";
    if(status == 0) {
      if (!ls->globals()["version"].is_nil() and !ls->globals()["version"]["cm_version"].is_nil()) {
        lua_cm_ver = ls->globals()["version"]["cm_version"].value<std::string>();
      }
    }

    // Show version information then quit if the -v option is used
    if(show_version) {
      std::cout << "Major version: " << MAJOR_VERSION << "\n" 
          "Code CM version: " << CM_VERSION << "\n";

      if (lua_cm_ver.length() > 0) {
        std::cout << "Lua config CM version: " << lua_cm_ver << "\n";
      }

      return 0;
    }

    boost::shared_ptr<L2FpConfigurationLua> config;
    if(Argc > 0) {
      std::string fname(Argv[0]);
      if(fname.size() > 4 && fname.substr(fname.size() - 4) == ".lua") {
        config.reset(new L2FpConfigurationLua(Argc, Argv));
        Logger::info() << "Level 2 Full Physics" << "\n"
                       << "Major Version: " << MAJOR_VERSION << "\n"
                       << "CM Version: " << CM_VERSION << "\n";
        if (lua_cm_ver.length() > 0) {
          Logger::info() << "Lua config CM version: " << lua_cm_ver << "\n";
        }
        Logger::info() << "Source directory: " << SOURCE_DIRECTORY << "\n"
                       << "Builder: " << BUILD_USER << "\n";
        char buf[PATH_MAX];
        Logger::info() << "Using Lua configuration file: "
                       << realpath(Argc > 0 ? Argv[0] :  "config.lua", buf)
                       << "\n";
      }
    }

    // Turn on signal for floating point errors. This helps catch more
    // obscure errors we have introduced in our code that might
    // introduce NAN or other floating point problems.
    //
    // Note:
    // 1. The Mac doesn't have this function, even though it is a C99
    //    function. We check for this during configuration.
    // 2. We Turn on FPE AFTER Lua because Lua puts numbers
    //    into integers first according to:
    //    http://lua-list.2524044.n2.nabble.com/FPE-with-using-gt-INT-MAX-numbers-in-Lua-td4985976.html
#ifdef HAVE_FEENABLEEXCEPT
    feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    struct sigaction act;
    act.sa_sigaction = floating_point_exception;
    sigemptyset(&act.sa_mask);
    act.sa_flags = SA_SIGINFO;    
    sigaction(SIGFPE, &act, (struct sigaction *)0);
    Logger::info() << "Turned on floating point exceptions\n";
#else
    Logger::info() << "This system does not support floating point exceptions\n";
#endif

    // Turn off the automatic ending of a program when GSL has an
    // error. Instead, we handle these as exceptions.
    gsl_set_error_handler_off();
    if(!config)
      throw Exception("Failed to load Lua configuration.");

    Logger::set_implementation(config->logger());
    time(&tm);
    Logger::info() << "Process started at: " << ctime(&tm) << "\n";
    // Set up output object
    boost::shared_ptr<Output> output;
    config->output(output, output_error);

    std::ostringstream tstream;
    tstream << "\n"
            << "==============================================\n"
            << "Setup:\n"
            << *config->forward_model() << "\n"
            << "==============================================\n";
    Logger::info() << tstream.str() << "\n";

    log_timing.write_to_log("Initialization");

    boost::shared_ptr<ConnorSolver> solver = config->solver();
    // Should have this merge with solver once we get this fully integrated
    boost::shared_ptr<IterativeSolver> iterative_solver = 
      config->iterative_solver();
    if (solver) {
      solver->add_observer(log_timing);

      boost::shared_ptr<InitialGuess> ig(config->initial_guess());
      blitz::Array<double, 1> initial_sv(ig->initial_guess());
      blitz::Array<double, 1> apriori_sv(ig->apriori());
      blitz::Array<double, 2> apriori_cov(ig->apriori_covariance());

      // Launch solver and iterate till an exit condition is reached
      bool have_sol = solver->solve(initial_sv, apriori_sv, apriori_cov);

      // Statevector isn't generally set to the final solution by
      // the Solver, so set it now.
      config->forward_model()->state_vector()->update_state(solver->x_solution(), 
                                                            solver->aposteriori_covariance());
      
      // Write output file
      output->write();

      if(have_sol)
        Logger::info() << "Found solution\n";
      else
        Logger::info() << "Failed to find solution\n";

      if(save_test) {        // Backdoor for generating test data.
        std::ofstream save(save_test_file.c_str());
        save <<
          "# This test data was generated by capturing the state of the ConnorSolver\n"
          "# after a run of l2_fp with the canned data used by \"make l2_fp_run\". You\n"
          "# can regenerate this by running l2_fp with the \"-t\" option.\n\n";
        save.precision(24);
        solver->to_stream(save);
      }
    } else if(iterative_solver) {
      // Note we don't have the setting of apriori and initial guess
      // like with solver. We want to get that in place, but for now
      // this gets set in the lua code and can't be changed here.
      Logger::info() << "SOLVER: A solver from the new solver class hierarchy\n"
;
      iterative_solver->solve();

      for( int i=0; i<=iterative_solver->num_accepted_steps(); i++ )
       Logger::info()
          << "   ========================================\n"
          << "   At point[" << i <<"] " 
	  << iterative_solver->accepted_points()[i] << "\n"
          << "   Cost[" << i << "] = " 
	  << iterative_solver->cost_at_accepted_points()[i] << "\n\n";

      boost::shared_ptr<MaxAPosteriori> map = config->max_a_posteriori();

      // Statevector isn't generally set to the final solution by
      // the Solver, so set it now.
      config->forward_model()->state_vector()->update_state
	(map->parameters(), map->a_posteriori_covariance());
      // Write output file
      output->write();      
      
      if(iterative_solver->status() == IterativeSolver::SUCCESS) 
	Logger::info() << "Found solution\n";
      else
	Logger::info() << "Failed to find solution\n";
    } else {
      // Run in forward model only mode if solver does not exist
      // Only calculate jacobians if they are to be written to 
      // the output file
      Logger::info() << "Running forward model only, no retrieval.\n";
      LuabindObject lua_config = config->lua_state().globals()["config"];
      bool jac_calc = true;
      if (not lua_config["write_jacobian"].is_nil())
        jac_calc = lua_config["write_jacobian"].value<bool>();
      Spectrum fm_radiance = config->forward_model()->radiance_all(jac_calc);

      // Store the radiance we just calculated into the output product
      output->register_data_source("/SpectralParameters/modeled_radiance", 
                                   fm_radiance.spectral_range().data());
      output->write();
    }

    log_timing.write_to_log("Final");
    time(&tm);
    Logger::info() << "Process ended at: " << ctime(&tm) << "\n";
    Logger::info() << "Bye bye\n";
  } catch(const Exception& e) {
    std::cerr << "Exception thrown by Full Physics code:\n"
              << e.what() << "\n"
              << "Back trace:\n" << boost::trace(e) << "\n";
    if(output_error)
      output_error->write_best_attempt();
    time(&tm);
    Logger::info() << "Process ended at: " << ctime(&tm) << "\n";
    log_timing.write_to_log("Final Error");
    return 1;
  } catch(const std::exception& e) {
    std::cerr << "System Exception thrown by Full Physics code:\n"
              << e.what() << "\n";
    if(output_error)
      output_error->write_best_attempt();
    time(&tm);
    Logger::info() << "Process ended at: " << ctime(&tm) << "\n";
    log_timing.write_to_log("Final Error");
    return 2;
  } catch(const H5::Exception& e) {
    std::cerr << "HDF 5 Exception thrown by Full Physics code:\n"
              << e.getDetailMsg() << "\n";
    if(output_error)
      output_error->write_best_attempt();
    time(&tm);
    Logger::info() << "Process ended at: " << ctime(&tm) << "\n";
    log_timing.write_to_log("Final Error");
    return 3;
  } catch(...) {
    std::cerr << "Unknown exception thrown\n";
    if(output_error)
      output_error->write_best_attempt();
    time(&tm);
    Logger::info() << "Process ended at: " << ctime(&tm) << "\n";
    log_timing.write_to_log("Final Error");
    return 4;
  }
  return 0;
}

