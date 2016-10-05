// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// This file was auto-generated

%include "common.i"
%{
#include "lidort_interface_masters.h"
%}
%import "lidort_interface_types.i"

namespace FullPhysics {



class Brdf_Linsup_Masters {

public:
  Brdf_Linsup_Masters(const int& thread_in);
  %python_attribute(brdf_sup_in, Brdf_Sup_Inputs&)
  %python_attribute(brdf_sup_inputstatus, Brdf_Input_Exception_Handling&)
  %python_attribute(thread, int)
  %python_attribute(brdf_sup_out, Brdf_Sup_Outputs&)
  %python_attribute(brdf_linsup_out, Brdf_Linsup_Outputs&)
  
  void read_config(const std::string& filnam_in);
  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in);
};


class Brdf_Sup_Masters {

public:
  Brdf_Sup_Masters(const int& thread_in);
  %python_attribute(brdf_sup_in, Brdf_Sup_Inputs&)
  %python_attribute(brdf_sup_inputstatus,Brdf_Input_Exception_Handling&)
  %python_attribute(thread,int)
  %python_attribute(brdf_sup_out,Brdf_Sup_Outputs&)
  
  void read_config(const std::string& filnam_in);
  void run(const bool& do_debug_restoration_in, const int& nmoments_input_in);
};


class Lidort_Lcs_Masters {

public:
  Lidort_Lcs_Masters(const int& thread_in);

  %python_attribute(thread,int)
  %python_attribute(lidort_fixin,Lidort_Fixed_Inputs&)
  %python_attribute(lidort_modin,Lidort_Modified_Inputs&)
  %python_attribute(lidort_sup,Lidort_Sup_Inout&)
  %python_attribute(lidort_out,Lidort_Outputs&)
  %python_attribute(lidort_linfixin,Lidort_Fixed_Lininputs&)
  %python_attribute(lidort_linmodin,Lidort_Modified_Lininputs&)
  %python_attribute(lidort_linsup,Lidort_Linsup_Inout&)
  %python_attribute(lidort_linout,Lidort_Linoutputs&)
  
  void run();
};


class Lidort_Lps_Masters {

public:
  Lidort_Lps_Masters(const int& thread_in);
  
  %python_attribute(thread,int)
  %python_attribute(lidort_fixin,Lidort_Fixed_Inputs&)
  %python_attribute(lidort_modin,Lidort_Modified_Inputs&)
  %python_attribute(lidort_sup,Lidort_Sup_Inout&)
  %python_attribute(lidort_out,Lidort_Outputs&)
  %python_attribute(lidort_linfixin,Lidort_Fixed_Lininputs&)
  %python_attribute(lidort_linmodin,Lidort_Modified_Lininputs&)
  %python_attribute(lidort_linsup,Lidort_Linsup_Inout&)
  %python_attribute(lidort_linout,Lidort_Linoutputs&)
  
  void run();
};


class Lidort_Inputs {

public:
  Lidort_Inputs();
  
  %python_attribute(lidort_fixin,Lidort_Fixed_Inputs&)
  %python_attribute(lidort_modin,Lidort_Modified_Inputs&)
  %python_attribute(lidort_inputstatus,Lidort_Input_Exception_Handling&)
  
  void read_config(const std::string& filnam_in);
};


class Lidort_Masters {

public:
  Lidort_Masters(const int& thread_in);
  
  %python_attribute(thread,int)
  %python_attribute(lidort_fixin,Lidort_Fixed_Inputs&)
  %python_attribute(lidort_modin,Lidort_Modified_Inputs&)
  %python_attribute(lidort_sup,Lidort_Sup_Inout&)
  %python_attribute(lidort_out,Lidort_Outputs&)
  
  void run();
};


class Lidort_Sup_Accessories {

public:
  Lidort_Sup_Accessories(boost::shared_ptr<Brdf_Sup_Inputs>& brdf_sup_in_in, boost::shared_ptr<Lidort_Fixed_Inputs>& lidort_fixin_in, boost::shared_ptr<Lidort_Modified_Inputs>& lidort_modin_in);
  
  %python_attribute(brdf_sup_in,Brdf_Sup_Inputs&)
  %python_attribute(lidort_fixin,Lidort_Fixed_Inputs&)
  %python_attribute(lidort_modin,Lidort_Modified_Inputs&)
  %python_attribute(lidort_brdfcheck_status,Lidort_Exception_Handling&)
  
  void brdf_input_checker();
};

}
