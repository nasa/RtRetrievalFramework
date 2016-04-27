#ifndef L2_FP_CONFIGURATION_LUA_H
#define L2_FP_CONFIGURATION_LUA_H
#include "l2_fp_configuration.h"
#include "lua_state.h"

namespace FullPhysics {
/****************************************************************//**
   This is an implementation of L2FpConfiguration that uses Lua. Note 
   that this class is very specific to the L2 Full Physics main, with 
   a minimal interface. A more general class LuaFile may be of more 
   interest to you unless you are specifically working with L2 main.

   The interface is purposely minimal, there are a handleful of global
   variables that the Lua file need to define. It can do this in any
   way it wishes. 

   In practice, our runs tend to be similiar. There is a file
   config_common.lua and base_config.lua that sets up a "standard"
   run, you can then create a configuration file that uses these two
   other files and just specifies what is different from this. But
   this is merely meant for convience, there is no requirement at all
   that you use these files. You can do anything you like as long as
   in the end you produce the set of global variables.

   The require variables are:

   \li logger - A LogImp
   \li forward_model - A ForwardModel
   \li solver - A ConnorSolver
   \li initial_guess - A InitialGuess
   \li number_pressure_level - Integer giving the number of pressure 
       levels. This is used to size the output file, so it should
       be the maximim pressure
   \li number_aerosol - Integer giving number of Aerosol particles. 
       Can be 0. This is used to size the output file.
   \li iteration_output - Boolean. True if we should write out iteration
       output
   \li register_output - A
       std::vector<boost::shared_ptr<RegisterOutputBase> > giving the
       list of output that should be generated. This list can empty if
       no output is desired. The Lua type for this is called 
       VectorRegisterOutput (since Lua doesn't have templates).
*******************************************************************/

class L2FpConfigurationLua : public L2FpConfiguration {
public:
//-----------------------------------------------------------------------
/// Read the given Lua configuration file. See the class description for
/// what this file needs to define.
//-----------------------------------------------------------------------

  L2FpConfigurationLua(const std::string& Fname, 
		       const std::string& Out_file = "out.h5")
    : ls(LuaState::load_file(Fname)), output_name_(Out_file) {}

//-----------------------------------------------------------------------
/// Take a LuaState that we've already generated the configuration
/// stuff (e.g., this was created in python script). See the class
/// description for what this state needs to define.
//-----------------------------------------------------------------------

  L2FpConfigurationLua(const boost::shared_ptr<LuaState>& Ls, 
		       const std::string& Out_file = "out.h5")
    : ls(Ls), output_name_(Out_file) {}
  
//-----------------------------------------------------------------------
/// Parse the arguments passed to the executable to set up
/// configuration. 
//-----------------------------------------------------------------------

  L2FpConfigurationLua(int Argc, char** Argv)
    : ls(LuaState::load_file(Argc > 0 ? Argv[0] :  "config.lua")), 
      output_name_(Argc > 1 ? Argv[1] : "out.h5")
  {
  }
  virtual ~L2FpConfigurationLua() {}
  virtual boost::shared_ptr<LogImp> logger() const
  {
    return ls->globals()["logger"].value_ptr<LogImp>();
  }
  virtual boost::shared_ptr<ForwardModel> forward_model() const 
  {
    return ls->globals()["forward_model"].value_ptr<ForwardModel>();
  }
  virtual boost::shared_ptr<ConnorSolver> solver() const
  {
    if (not ls->globals()["solver"].is_nil())
      return ls->globals()["solver"].value_ptr<ConnorSolver>();
    else
      return boost::shared_ptr<ConnorSolver>();
  }
  virtual boost::shared_ptr<InitialGuess> initial_guess() const
  {
    return ls->globals()["initial_guess"].value_ptr<InitialGuess>();
  }
  virtual void output(boost::shared_ptr<Output>& Regular_output,
		      boost::shared_ptr<Output>& Error_output) const;
  virtual void print(std::ostream& Os) const { Os << "L2FpConfigurationLua"; }
  LuaState& lua_state() {return *ls;}

  virtual boost::shared_ptr<IterativeSolver> iterative_solver() const
  {
    if (not ls->globals()["iterative_solver"].is_nil())
      return ls->globals()["iterative_solver"].value_ptr<IterativeSolver>();
    return boost::shared_ptr<IterativeSolver>();
  }

  virtual boost::shared_ptr<MaxAPosteriori> max_a_posteriori() const
  {
    if(not ls->globals()["stat_method_map"].is_nil())
      return ls->globals()["stat_method_map"].value_ptr<MaxAPosteriori>();
    return boost::shared_ptr<MaxAPosteriori>();
  }

//-----------------------------------------------------------------------
/// Output file name.
//-----------------------------------------------------------------------

  const std::string& output_name() const {return output_name_;}
  void output_name(const std::string& F) { output_name_ = F; }
private:
  boost::shared_ptr<LuaState> ls;
  std::string output_name_;
  mutable boost::shared_ptr<Output> out_iteration;
};
}
#endif
