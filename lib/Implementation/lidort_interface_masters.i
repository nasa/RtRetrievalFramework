// This file was auto-generated

%include "common.i"

%{
#include "lidort_interface_masters.h"
%}

%import "lidort_interface_types.i"

%fp_shared_ptr(FullPhysics::Brdf_Lin_Sup_Masters);
%fp_shared_ptr(FullPhysics::Brdf_Sup_Masters);
%fp_shared_ptr(FullPhysics::Lidort_Inputs);
%fp_shared_ptr(FullPhysics::Lidort_Masters);
%fp_shared_ptr(FullPhysics::Lidort_L_Inputs);
%fp_shared_ptr(FullPhysics::Lidort_Lcs_Masters);
%fp_shared_ptr(FullPhysics::Lidort_Lps_Masters);
%fp_shared_ptr(FullPhysics::Lidort_Brdf_Sup_Accessories);

namespace FullPhysics {



class Brdf_Lin_Sup_Masters {

public:
  Brdf_Lin_Sup_Masters();
  virtual ~Brdf_Lin_Sup_Masters();
  std::string print_to_string() const;

  %python_attribute(brdf_sup_in, Brdf_Sup_Inputs&)
  %python_attribute(brdf_linsup_in, Brdf_Linsup_Inputs&)
  %python_attribute(brdf_sup_inputstatus, Brdf_Input_Exception_Handling&)
  %python_attribute(brdf_sup_out, Brdf_Sup_Outputs&)
  %python_attribute(brdf_linsup_out, Brdf_Linsup_Outputs&)
  %python_attribute(brdf_sup_outputstatus, Brdf_Output_Exception_Handling&)
  
  void read_config(const std::string& filnam_in);
  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in);
};


class Brdf_Sup_Masters {

public:
  Brdf_Sup_Masters();
  virtual ~Brdf_Sup_Masters();
  std::string print_to_string() const;

  %python_attribute(brdf_sup_in, Brdf_Sup_Inputs&)
  %python_attribute(brdf_sup_inputstatus, Brdf_Input_Exception_Handling&)
  %python_attribute(brdf_sup_out, Brdf_Sup_Outputs&)
  %python_attribute(brdf_sup_outputstatus, Brdf_Output_Exception_Handling&)
  
  void read_config(const std::string& filnam_in);
  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in);
};


class Lidort_Inputs {

public:
  Lidort_Inputs();
  virtual ~Lidort_Inputs();
  std::string print_to_string() const;

  %python_attribute(lidort_sup, Lidort_Sup_Inout&)
  %python_attribute(lidort_fixin, Lidort_Fixed_Inputs&)
  %python_attribute(lidort_modin, Lidort_Modified_Inputs&)
  %python_attribute(lidort_inputstatus, Lidort_Input_Exception_Handling&)
  
  void brdf_sup_init();
  void read_config(const std::string& filnam_in);
  void sleave_sup_init();
  void ss_sup_init();
  void sup_init();
};


class Lidort_Masters {

public:
  Lidort_Masters();
  virtual ~Lidort_Masters();
  std::string print_to_string() const;

  %python_attribute(lidort_fixin, Lidort_Fixed_Inputs&)
  %python_attribute(lidort_modin, Lidort_Modified_Inputs&)
  %python_attribute(lidort_sup, Lidort_Sup_Inout&)
  %python_attribute(lidort_out, Lidort_Outputs&)
  
  void run(const bool& do_debug_input_in);
};


class Lidort_L_Inputs {

public:
  Lidort_L_Inputs();
  virtual ~Lidort_L_Inputs();
  std::string print_to_string() const;

  %python_attribute(lidort_linsup, Lidort_Linsup_Inout&)
  %python_attribute(lidort_fixin, Lidort_Fixed_Inputs&)
  %python_attribute(lidort_modin, Lidort_Modified_Inputs&)
  %python_attribute(lidort_linfixin, Lidort_Fixed_Lininputs&)
  %python_attribute(lidort_linmodin, Lidort_Modified_Lininputs&)
  %python_attribute(lidort_inputstatus, Lidort_Input_Exception_Handling&)
  
  void brdf_linsup_init();
  void read_config(const std::string& filnam_in);
  void linsup_init();
  void sleave_linsup_init();
  void ss_linsup_init();
};


class Lidort_Lcs_Masters {

public:
  Lidort_Lcs_Masters();
  virtual ~Lidort_Lcs_Masters();
  std::string print_to_string() const;

  %python_attribute(lidort_fixin, Lidort_Fixed_Inputs&)
  %python_attribute(lidort_modin, Lidort_Modified_Inputs&)
  %python_attribute(lidort_sup, Lidort_Sup_Inout&)
  %python_attribute(lidort_out, Lidort_Outputs&)
  %python_attribute(lidort_linfixin, Lidort_Fixed_Lininputs&)
  %python_attribute(lidort_linmodin, Lidort_Modified_Lininputs&)
  %python_attribute(lidort_linsup, Lidort_Linsup_Inout&)
  %python_attribute(lidort_linout, Lidort_Linoutputs&)
  
  void run(const bool& do_debug_input_in);
};


class Lidort_Lps_Masters {

public:
  Lidort_Lps_Masters();
  virtual ~Lidort_Lps_Masters();
  std::string print_to_string() const;

  %python_attribute(lidort_fixin, Lidort_Fixed_Inputs&)
  %python_attribute(lidort_modin, Lidort_Modified_Inputs&)
  %python_attribute(lidort_sup, Lidort_Sup_Inout&)
  %python_attribute(lidort_out, Lidort_Outputs&)
  %python_attribute(lidort_linfixin, Lidort_Fixed_Lininputs&)
  %python_attribute(lidort_linmodin, Lidort_Modified_Lininputs&)
  %python_attribute(lidort_linsup, Lidort_Linsup_Inout&)
  %python_attribute(lidort_linout, Lidort_Linoutputs&)
  
  void run(const bool& do_debug_input_in);
};


class Lidort_Brdf_Sup_Accessories {

public:
  Lidort_Brdf_Sup_Accessories(boost::shared_ptr<Brdf_Sup_Inputs>& brdf_sup_in_in, boost::shared_ptr<Lidort_Fixed_Inputs>& lidort_fixin_in, boost::shared_ptr<Lidort_Modified_Inputs>& lidort_modin_in);
  virtual ~Lidort_Brdf_Sup_Accessories();
  std::string print_to_string() const;

  %python_attribute(brdf_sup_in, Brdf_Sup_Inputs&)
  %python_attribute(lidort_fixin, Lidort_Fixed_Inputs&)
  %python_attribute(lidort_modin, Lidort_Modified_Inputs&)
  %python_attribute(lidort_brdfcheck_status, Lidort_Exception_Handling&)
  %python_attribute(brdf_sup_out, Brdf_Sup_Outputs&)
  %python_attribute(lidort_sup, Lidort_Sup_Inout&)
  
  void brdf_input_check();
  void brdf_input_check_error(const std::string& errorfile_in);
  void set_brdf_inputs();
};

}