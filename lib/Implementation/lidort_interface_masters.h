#ifndef LIDORT_INTERFACE_MASTERS_H
#define LIDORT_INTERFACE_MASTERS_H

#include <iostream>
#include <blitz/array.h>

#include "fp_exception.h"
#include "lidort_interface_types.h"

/* This file was auto-generated */

namespace FullPhysics {

//-----------------------------------------------------------------------
// Links to module: "brdf_linsup_masters_m" in file: "brdf_lin_sup_masters.F90"
//-----------------------------------------------------------------------

extern "C" {
  void brdf_linsup_masters_m_read_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** brdf_linsup_in_in, void** brdf_sup_inputstatus_in, void** brdf_sup_out_in, void** brdf_linsup_out_in);
  void brdf_linsup_masters_m_write_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** brdf_linsup_in_in, void** brdf_sup_inputstatus_in, void** brdf_sup_out_in, void** brdf_linsup_out_in);
  void brdf_linsup_masters_m_brdf_lin_inputmaster_wrap(const int* filnam_in_len, const char* filnam_in, void** brdf_sup_in_in, void** brdf_linsup_in_in, void** brdf_sup_inputstatus_in);
  void brdf_linsup_masters_m_brdf_lin_mainmaster_wrap(const int* thread_in, const bool* do_debug_restoration_in, const int* nmoments_input_in, void** brdf_sup_in_in, void** brdf_linsup_in_in, void** brdf_sup_out_in, void** brdf_linsup_out_in);
}

class Brdf_Linsup_Masters {

public:
  Brdf_Linsup_Masters(const int& thread_in) : thread_(thread_in) 
  { 
    
    // Initialize type pointers
    brdf_sup_in_.reset( new Brdf_Sup_Inputs() );
    brdf_linsup_in_.reset( new Brdf_Linsup_Inputs() );
    brdf_sup_inputstatus_.reset( new Brdf_Input_Exception_Handling() );
    brdf_sup_out_.reset( new Brdf_Sup_Outputs() );
    brdf_linsup_out_.reset( new Brdf_Linsup_Outputs() );
    
  }

  Brdf_Sup_Inputs& brdf_sup_in() {
    return *brdf_sup_in_;
  }
  const Brdf_Sup_Inputs& brdf_sup_in() const {
    return *brdf_sup_in_;
  }
  boost::shared_ptr<Brdf_Sup_Inputs>& brdf_sup_in_ptr() {
    return brdf_sup_in_;
  }
  

  Brdf_Linsup_Inputs& brdf_linsup_in() {
    return *brdf_linsup_in_;
  }
  const Brdf_Linsup_Inputs& brdf_linsup_in() const {
    return *brdf_linsup_in_;
  }
  boost::shared_ptr<Brdf_Linsup_Inputs>& brdf_linsup_in_ptr() {
    return brdf_linsup_in_;
  }
  

  Brdf_Input_Exception_Handling& brdf_sup_inputstatus() {
    return *brdf_sup_inputstatus_;
  }
  const Brdf_Input_Exception_Handling& brdf_sup_inputstatus() const {
    return *brdf_sup_inputstatus_;
  }
  boost::shared_ptr<Brdf_Input_Exception_Handling>& brdf_sup_inputstatus_ptr() {
    return brdf_sup_inputstatus_;
  }
  

  const int& thread() const {
    return thread_;
  }
  

  Brdf_Sup_Outputs& brdf_sup_out() {
    return *brdf_sup_out_;
  }
  const Brdf_Sup_Outputs& brdf_sup_out() const {
    return *brdf_sup_out_;
  }
  boost::shared_ptr<Brdf_Sup_Outputs>& brdf_sup_out_ptr() {
    return brdf_sup_out_;
  }
  

  Brdf_Linsup_Outputs& brdf_linsup_out() {
    return *brdf_linsup_out_;
  }
  const Brdf_Linsup_Outputs& brdf_linsup_out() const {
    return *brdf_linsup_out_;
  }
  boost::shared_ptr<Brdf_Linsup_Outputs>& brdf_linsup_out_ptr() {
    return brdf_linsup_out_;
  }
  

  
  void read_config(const std::string& filnam_in) {
    const char* filnam_lcl = filnam_in.c_str();
    int filnam_in_len = (int) filnam_in.size();
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_linsup_in_lcl = brdf_linsup_in_->fortran_type_ptr();
    void* brdf_sup_inputstatus_lcl = brdf_sup_inputstatus_->fortran_type_ptr();
    
    brdf_linsup_masters_m_brdf_lin_inputmaster_wrap(&filnam_in_len, filnam_lcl, &brdf_sup_in_lcl, &brdf_linsup_in_lcl, &brdf_sup_inputstatus_lcl);
    

  }
void run(const bool& do_debug_restoration_in, const int& nmoments_input_in) {
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_linsup_in_lcl = brdf_linsup_in_->fortran_type_ptr();
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();
    void* brdf_linsup_out_lcl = brdf_linsup_out_->fortran_type_ptr();
    
    brdf_linsup_masters_m_brdf_lin_mainmaster_wrap(&thread_, &do_debug_restoration_in, &nmoments_input_in, &brdf_sup_in_lcl, &brdf_linsup_in_lcl, &brdf_sup_out_lcl, &brdf_linsup_out_lcl);
    

  }
 
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_linsup_in_lcl = brdf_linsup_in_->fortran_type_ptr();
    void* brdf_sup_inputstatus_lcl = brdf_sup_inputstatus_->fortran_type_ptr();
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();
    void* brdf_linsup_out_lcl = brdf_linsup_out_->fortran_type_ptr();

    brdf_linsup_masters_m_read_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &brdf_linsup_in_lcl, &brdf_sup_inputstatus_lcl, &brdf_sup_out_lcl, &brdf_linsup_out_lcl);
    
  }
   
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_linsup_in_lcl = brdf_linsup_in_->fortran_type_ptr();
    void* brdf_sup_inputstatus_lcl = brdf_sup_inputstatus_->fortran_type_ptr();
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();
    void* brdf_linsup_out_lcl = brdf_linsup_out_->fortran_type_ptr();

    brdf_linsup_masters_m_write_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &brdf_linsup_in_lcl, &brdf_sup_inputstatus_lcl, &brdf_sup_out_lcl, &brdf_linsup_out_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Brdf_Linsup_Masters &obj) {
    output_stream << "Brdf_Linsup_Masters:" << std::endl
      << "         brdf_sup_in: " << obj.brdf_sup_in()  << std::endl
      << "      brdf_linsup_in: " << obj.brdf_linsup_in()  << std::endl
      << "brdf_sup_inputstatus: " << obj.brdf_sup_inputstatus()  << std::endl
      << "              thread: " << obj.thread()  << std::endl
      << "        brdf_sup_out: " << obj.brdf_sup_out()  << std::endl
      << "     brdf_linsup_out: " << obj.brdf_linsup_out()  << std::endl;
    return output_stream;
  }

private:
  boost::shared_ptr<Brdf_Sup_Inputs> brdf_sup_in_;
  boost::shared_ptr<Brdf_Linsup_Inputs> brdf_linsup_in_;
  boost::shared_ptr<Brdf_Input_Exception_Handling> brdf_sup_inputstatus_;
  int thread_;
  boost::shared_ptr<Brdf_Sup_Outputs> brdf_sup_out_;
  boost::shared_ptr<Brdf_Linsup_Outputs> brdf_linsup_out_;
};

//-----------------------------------------------------------------------
// Links to module: "brdf_sup_masters_m" in file: "brdf_sup_masters.F90"
//-----------------------------------------------------------------------

extern "C" {
  void brdf_sup_masters_m_read_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** brdf_sup_inputstatus_in, void** brdf_sup_out_in);
  void brdf_sup_masters_m_write_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** brdf_sup_inputstatus_in, void** brdf_sup_out_in);
  void brdf_sup_masters_m_brdf_inputmaster_wrap(const int* filnam_in_len, const char* filnam_in, void** brdf_sup_in_in, void** brdf_sup_inputstatus_in);
  void brdf_sup_masters_m_brdf_mainmaster_wrap(const int* thread_in, const bool* do_debug_restoration_in, const int* nmoments_input_in, void** brdf_sup_in_in, void** brdf_sup_out_in);
}

class Brdf_Sup_Masters {

public:
  Brdf_Sup_Masters(const int& thread_in) : thread_(thread_in) 
  { 
    
    // Initialize type pointers
    brdf_sup_in_.reset( new Brdf_Sup_Inputs() );
    brdf_sup_inputstatus_.reset( new Brdf_Input_Exception_Handling() );
    brdf_sup_out_.reset( new Brdf_Sup_Outputs() );
    
  }

  Brdf_Sup_Inputs& brdf_sup_in() {
    return *brdf_sup_in_;
  }
  const Brdf_Sup_Inputs& brdf_sup_in() const {
    return *brdf_sup_in_;
  }
  boost::shared_ptr<Brdf_Sup_Inputs>& brdf_sup_in_ptr() {
    return brdf_sup_in_;
  }
  

  Brdf_Input_Exception_Handling& brdf_sup_inputstatus() {
    return *brdf_sup_inputstatus_;
  }
  const Brdf_Input_Exception_Handling& brdf_sup_inputstatus() const {
    return *brdf_sup_inputstatus_;
  }
  boost::shared_ptr<Brdf_Input_Exception_Handling>& brdf_sup_inputstatus_ptr() {
    return brdf_sup_inputstatus_;
  }
  

  const int& thread() const {
    return thread_;
  }
  

  Brdf_Sup_Outputs& brdf_sup_out() {
    return *brdf_sup_out_;
  }
  const Brdf_Sup_Outputs& brdf_sup_out() const {
    return *brdf_sup_out_;
  }
  boost::shared_ptr<Brdf_Sup_Outputs>& brdf_sup_out_ptr() {
    return brdf_sup_out_;
  }
  

  
  void read_config(const std::string& filnam_in) {
    const char* filnam_lcl = filnam_in.c_str();
    int filnam_in_len = (int) filnam_in.size();
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_sup_inputstatus_lcl = brdf_sup_inputstatus_->fortran_type_ptr();
    
    brdf_sup_masters_m_brdf_inputmaster_wrap(&filnam_in_len, filnam_lcl, &brdf_sup_in_lcl, &brdf_sup_inputstatus_lcl);
    

  }
void run(const bool& do_debug_restoration_in, const int& nmoments_input_in) {
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();
    
    brdf_sup_masters_m_brdf_mainmaster_wrap(&thread_, &do_debug_restoration_in, &nmoments_input_in, &brdf_sup_in_lcl, &brdf_sup_out_lcl);
    

  }
 
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_sup_inputstatus_lcl = brdf_sup_inputstatus_->fortran_type_ptr();
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();

    brdf_sup_masters_m_read_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &brdf_sup_inputstatus_lcl, &brdf_sup_out_lcl);
    
  }
   
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_sup_inputstatus_lcl = brdf_sup_inputstatus_->fortran_type_ptr();
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();

    brdf_sup_masters_m_write_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &brdf_sup_inputstatus_lcl, &brdf_sup_out_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Brdf_Sup_Masters &obj) {
    output_stream << "Brdf_Sup_Masters:" << std::endl
      << "         brdf_sup_in: " << obj.brdf_sup_in()  << std::endl
      << "brdf_sup_inputstatus: " << obj.brdf_sup_inputstatus()  << std::endl
      << "              thread: " << obj.thread()  << std::endl
      << "        brdf_sup_out: " << obj.brdf_sup_out()  << std::endl;
    return output_stream;
  }

private:
  boost::shared_ptr<Brdf_Sup_Inputs> brdf_sup_in_;
  boost::shared_ptr<Brdf_Input_Exception_Handling> brdf_sup_inputstatus_;
  int thread_;
  boost::shared_ptr<Brdf_Sup_Outputs> brdf_sup_out_;
};

//-----------------------------------------------------------------------
// Links to module: "lidort_lcs_masters" in file: "lidort_lcs_masters.F90"
//-----------------------------------------------------------------------

extern "C" {
  void lcs_masters_read_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
  void lcs_masters_write_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
  void lcs_masters_lcs_master_wrap(const int* thread_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
}

class Lidort_Lcs_Masters {

public:
  Lidort_Lcs_Masters(const int& thread_in) : thread_(thread_in) 
  { 
    
    // Initialize type pointers
    lidort_fixin_.reset( new Lidort_Fixed_Inputs() );
    lidort_modin_.reset( new Lidort_Modified_Inputs() );
    lidort_sup_.reset( new Lidort_Sup_Inout() );
    lidort_out_.reset( new Lidort_Outputs() );
    lidort_linfixin_.reset( new Lidort_Fixed_Lininputs() );
    lidort_linmodin_.reset( new Lidort_Modified_Lininputs() );
    lidort_linsup_.reset( new Lidort_Linsup_Inout() );
    lidort_linout_.reset( new Lidort_Linoutputs() );
    
  }

  const int& thread() const {
    return thread_;
  }
  

  Lidort_Fixed_Inputs& lidort_fixin() {
    return *lidort_fixin_;
  }
  const Lidort_Fixed_Inputs& lidort_fixin() const {
    return *lidort_fixin_;
  }
  boost::shared_ptr<Lidort_Fixed_Inputs>& lidort_fixin_ptr() {
    return lidort_fixin_;
  }
  void lidort_fixin(Lidort_Fixed_Inputs& lidort_fixin_in) {
    void* src_ptr = lidort_fixin_in.fortran_type_ptr();
    void* dst_ptr = lidort_fixin_->fortran_type_ptr();
    lidort_fixed_inputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Modified_Inputs& lidort_modin() {
    return *lidort_modin_;
  }
  const Lidort_Modified_Inputs& lidort_modin() const {
    return *lidort_modin_;
  }
  boost::shared_ptr<Lidort_Modified_Inputs>& lidort_modin_ptr() {
    return lidort_modin_;
  }
  void lidort_modin(Lidort_Modified_Inputs& lidort_modin_in) {
    void* src_ptr = lidort_modin_in.fortran_type_ptr();
    void* dst_ptr = lidort_modin_->fortran_type_ptr();
    lidort_modified_inputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Sup_Inout& lidort_sup() {
    return *lidort_sup_;
  }
  const Lidort_Sup_Inout& lidort_sup() const {
    return *lidort_sup_;
  }
  boost::shared_ptr<Lidort_Sup_Inout>& lidort_sup_ptr() {
    return lidort_sup_;
  }
  void lidort_sup(Lidort_Sup_Inout& lidort_sup_in) {
    void* src_ptr = lidort_sup_in.fortran_type_ptr();
    void* dst_ptr = lidort_sup_->fortran_type_ptr();
    lidort_sup_inout_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Outputs& lidort_out() {
    return *lidort_out_;
  }
  const Lidort_Outputs& lidort_out() const {
    return *lidort_out_;
  }
  boost::shared_ptr<Lidort_Outputs>& lidort_out_ptr() {
    return lidort_out_;
  }
  void lidort_out(Lidort_Outputs& lidort_out_in) {
    void* src_ptr = lidort_out_in.fortran_type_ptr();
    void* dst_ptr = lidort_out_->fortran_type_ptr();
    lidort_outputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Fixed_Lininputs& lidort_linfixin() {
    return *lidort_linfixin_;
  }
  const Lidort_Fixed_Lininputs& lidort_linfixin() const {
    return *lidort_linfixin_;
  }
  boost::shared_ptr<Lidort_Fixed_Lininputs>& lidort_linfixin_ptr() {
    return lidort_linfixin_;
  }
  void lidort_linfixin(Lidort_Fixed_Lininputs& lidort_linfixin_in) {
    void* src_ptr = lidort_linfixin_in.fortran_type_ptr();
    void* dst_ptr = lidort_linfixin_->fortran_type_ptr();
    lidort_fixed_lininputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Modified_Lininputs& lidort_linmodin() {
    return *lidort_linmodin_;
  }
  const Lidort_Modified_Lininputs& lidort_linmodin() const {
    return *lidort_linmodin_;
  }
  boost::shared_ptr<Lidort_Modified_Lininputs>& lidort_linmodin_ptr() {
    return lidort_linmodin_;
  }
  void lidort_linmodin(Lidort_Modified_Lininputs& lidort_linmodin_in) {
    void* src_ptr = lidort_linmodin_in.fortran_type_ptr();
    void* dst_ptr = lidort_linmodin_->fortran_type_ptr();
    lidort_modified_lininputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Linsup_Inout& lidort_linsup() {
    return *lidort_linsup_;
  }
  const Lidort_Linsup_Inout& lidort_linsup() const {
    return *lidort_linsup_;
  }
  boost::shared_ptr<Lidort_Linsup_Inout>& lidort_linsup_ptr() {
    return lidort_linsup_;
  }
  void lidort_linsup(Lidort_Linsup_Inout& lidort_linsup_in) {
    void* src_ptr = lidort_linsup_in.fortran_type_ptr();
    void* dst_ptr = lidort_linsup_->fortran_type_ptr();
    lidort_linsup_inout_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Linoutputs& lidort_linout() {
    return *lidort_linout_;
  }
  const Lidort_Linoutputs& lidort_linout() const {
    return *lidort_linout_;
  }
  boost::shared_ptr<Lidort_Linoutputs>& lidort_linout_ptr() {
    return lidort_linout_;
  }
  void lidort_linout(Lidort_Linoutputs& lidort_linout_in) {
    void* src_ptr = lidort_linout_in.fortran_type_ptr();
    void* dst_ptr = lidort_linout_->fortran_type_ptr();
    lidort_linoutputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  
  void run() {
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();
    void* lidort_linfixin_lcl = lidort_linfixin_->fortran_type_ptr();
    void* lidort_linmodin_lcl = lidort_linmodin_->fortran_type_ptr();
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    void* lidort_linout_lcl = lidort_linout_->fortran_type_ptr();
    
    lcs_masters_lcs_master_wrap(&thread_, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
    Lidort_Pars lid_pars = Lidort_Pars::instance();
    if( lidort_out().status().ts_status_inputcheck() != lid_pars.lidort_success ||
        lidort_out().status().ts_status_calculation() != lid_pars.lidort_success ) {
       std::stringstream err_msg;
       err_msg << "LIDORT Error at " << __FILE__ << ":" << __LINE__ << std::endl;
       // Output the full details of the error message to stderr since the exception
       // class may truncate the message
       std::cerr << err_msg.str();
       std::cerr << lidort_out().status();
       throw Exception(err_msg.str());
    }

  }
 
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();
    void* lidort_linfixin_lcl = lidort_linfixin_->fortran_type_ptr();
    void* lidort_linmodin_lcl = lidort_linmodin_->fortran_type_ptr();
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    void* lidort_linout_lcl = lidort_linout_->fortran_type_ptr();

    lcs_masters_read_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
  }
   
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();
    void* lidort_linfixin_lcl = lidort_linfixin_->fortran_type_ptr();
    void* lidort_linmodin_lcl = lidort_linmodin_->fortran_type_ptr();
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    void* lidort_linout_lcl = lidort_linout_->fortran_type_ptr();

    lcs_masters_write_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_Lcs_Masters &obj) {
    output_stream << "Lidort_Lcs_Masters:" << std::endl
      << "         thread: " << obj.thread()  << std::endl
      << "   lidort_fixin: " << obj.lidort_fixin()  << std::endl
      << "   lidort_modin: " << obj.lidort_modin()  << std::endl
      << "     lidort_sup: " << obj.lidort_sup()  << std::endl
      << "     lidort_out: " << obj.lidort_out()  << std::endl
      << "lidort_linfixin: " << obj.lidort_linfixin()  << std::endl
      << "lidort_linmodin: " << obj.lidort_linmodin()  << std::endl
      << "  lidort_linsup: " << obj.lidort_linsup()  << std::endl
      << "  lidort_linout: " << obj.lidort_linout()  << std::endl;
    return output_stream;
  }

private:
  int thread_;
  boost::shared_ptr<Lidort_Fixed_Inputs> lidort_fixin_;
  boost::shared_ptr<Lidort_Modified_Inputs> lidort_modin_;
  boost::shared_ptr<Lidort_Sup_Inout> lidort_sup_;
  boost::shared_ptr<Lidort_Outputs> lidort_out_;
  boost::shared_ptr<Lidort_Fixed_Lininputs> lidort_linfixin_;
  boost::shared_ptr<Lidort_Modified_Lininputs> lidort_linmodin_;
  boost::shared_ptr<Lidort_Linsup_Inout> lidort_linsup_;
  boost::shared_ptr<Lidort_Linoutputs> lidort_linout_;
};

//-----------------------------------------------------------------------
// Links to module: "lidort_lps_masters" in file: "lidort_lps_masters.F90"
//-----------------------------------------------------------------------

extern "C" {
  void lps_masters_read_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
  void lps_masters_write_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
  void lps_masters_lps_master_wrap(const int* thread_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
}

class Lidort_Lps_Masters {

public:
  Lidort_Lps_Masters(const int& thread_in) : thread_(thread_in) 
  { 
    
    // Initialize type pointers
    lidort_fixin_.reset( new Lidort_Fixed_Inputs() );
    lidort_modin_.reset( new Lidort_Modified_Inputs() );
    lidort_sup_.reset( new Lidort_Sup_Inout() );
    lidort_out_.reset( new Lidort_Outputs() );
    lidort_linfixin_.reset( new Lidort_Fixed_Lininputs() );
    lidort_linmodin_.reset( new Lidort_Modified_Lininputs() );
    lidort_linsup_.reset( new Lidort_Linsup_Inout() );
    lidort_linout_.reset( new Lidort_Linoutputs() );
    
  }

  const int& thread() const {
    return thread_;
  }
  

  Lidort_Fixed_Inputs& lidort_fixin() {
    return *lidort_fixin_;
  }
  const Lidort_Fixed_Inputs& lidort_fixin() const {
    return *lidort_fixin_;
  }
  boost::shared_ptr<Lidort_Fixed_Inputs>& lidort_fixin_ptr() {
    return lidort_fixin_;
  }
  void lidort_fixin(Lidort_Fixed_Inputs& lidort_fixin_in) {
    void* src_ptr = lidort_fixin_in.fortran_type_ptr();
    void* dst_ptr = lidort_fixin_->fortran_type_ptr();
    lidort_fixed_inputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Modified_Inputs& lidort_modin() {
    return *lidort_modin_;
  }
  const Lidort_Modified_Inputs& lidort_modin() const {
    return *lidort_modin_;
  }
  boost::shared_ptr<Lidort_Modified_Inputs>& lidort_modin_ptr() {
    return lidort_modin_;
  }
  void lidort_modin(Lidort_Modified_Inputs& lidort_modin_in) {
    void* src_ptr = lidort_modin_in.fortran_type_ptr();
    void* dst_ptr = lidort_modin_->fortran_type_ptr();
    lidort_modified_inputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Sup_Inout& lidort_sup() {
    return *lidort_sup_;
  }
  const Lidort_Sup_Inout& lidort_sup() const {
    return *lidort_sup_;
  }
  boost::shared_ptr<Lidort_Sup_Inout>& lidort_sup_ptr() {
    return lidort_sup_;
  }
  void lidort_sup(Lidort_Sup_Inout& lidort_sup_in) {
    void* src_ptr = lidort_sup_in.fortran_type_ptr();
    void* dst_ptr = lidort_sup_->fortran_type_ptr();
    lidort_sup_inout_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Outputs& lidort_out() {
    return *lidort_out_;
  }
  const Lidort_Outputs& lidort_out() const {
    return *lidort_out_;
  }
  boost::shared_ptr<Lidort_Outputs>& lidort_out_ptr() {
    return lidort_out_;
  }
  void lidort_out(Lidort_Outputs& lidort_out_in) {
    void* src_ptr = lidort_out_in.fortran_type_ptr();
    void* dst_ptr = lidort_out_->fortran_type_ptr();
    lidort_outputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Fixed_Lininputs& lidort_linfixin() {
    return *lidort_linfixin_;
  }
  const Lidort_Fixed_Lininputs& lidort_linfixin() const {
    return *lidort_linfixin_;
  }
  boost::shared_ptr<Lidort_Fixed_Lininputs>& lidort_linfixin_ptr() {
    return lidort_linfixin_;
  }
  void lidort_linfixin(Lidort_Fixed_Lininputs& lidort_linfixin_in) {
    void* src_ptr = lidort_linfixin_in.fortran_type_ptr();
    void* dst_ptr = lidort_linfixin_->fortran_type_ptr();
    lidort_fixed_lininputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Modified_Lininputs& lidort_linmodin() {
    return *lidort_linmodin_;
  }
  const Lidort_Modified_Lininputs& lidort_linmodin() const {
    return *lidort_linmodin_;
  }
  boost::shared_ptr<Lidort_Modified_Lininputs>& lidort_linmodin_ptr() {
    return lidort_linmodin_;
  }
  void lidort_linmodin(Lidort_Modified_Lininputs& lidort_linmodin_in) {
    void* src_ptr = lidort_linmodin_in.fortran_type_ptr();
    void* dst_ptr = lidort_linmodin_->fortran_type_ptr();
    lidort_modified_lininputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Linsup_Inout& lidort_linsup() {
    return *lidort_linsup_;
  }
  const Lidort_Linsup_Inout& lidort_linsup() const {
    return *lidort_linsup_;
  }
  boost::shared_ptr<Lidort_Linsup_Inout>& lidort_linsup_ptr() {
    return lidort_linsup_;
  }
  void lidort_linsup(Lidort_Linsup_Inout& lidort_linsup_in) {
    void* src_ptr = lidort_linsup_in.fortran_type_ptr();
    void* dst_ptr = lidort_linsup_->fortran_type_ptr();
    lidort_linsup_inout_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Linoutputs& lidort_linout() {
    return *lidort_linout_;
  }
  const Lidort_Linoutputs& lidort_linout() const {
    return *lidort_linout_;
  }
  boost::shared_ptr<Lidort_Linoutputs>& lidort_linout_ptr() {
    return lidort_linout_;
  }
  void lidort_linout(Lidort_Linoutputs& lidort_linout_in) {
    void* src_ptr = lidort_linout_in.fortran_type_ptr();
    void* dst_ptr = lidort_linout_->fortran_type_ptr();
    lidort_linoutputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  
  void run() {
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();
    void* lidort_linfixin_lcl = lidort_linfixin_->fortran_type_ptr();
    void* lidort_linmodin_lcl = lidort_linmodin_->fortran_type_ptr();
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    void* lidort_linout_lcl = lidort_linout_->fortran_type_ptr();
    
    lps_masters_lps_master_wrap(&thread_, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
    Lidort_Pars lid_pars = Lidort_Pars::instance();
    if( lidort_out().status().ts_status_inputcheck() != lid_pars.lidort_success ||
        lidort_out().status().ts_status_calculation() != lid_pars.lidort_success ) {
       std::stringstream err_msg;
       err_msg << "LIDORT Error at " << __FILE__ << ":" << __LINE__ << std::endl;
       // Output the full details of the error message to stderr since the exception
       // class may truncate the message
       std::cerr << err_msg.str();
       std::cerr << lidort_out().status();
       throw Exception(err_msg.str());
    }

  }
 
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();
    void* lidort_linfixin_lcl = lidort_linfixin_->fortran_type_ptr();
    void* lidort_linmodin_lcl = lidort_linmodin_->fortran_type_ptr();
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    void* lidort_linout_lcl = lidort_linout_->fortran_type_ptr();

    lps_masters_read_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
  }
   
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();
    void* lidort_linfixin_lcl = lidort_linfixin_->fortran_type_ptr();
    void* lidort_linmodin_lcl = lidort_linmodin_->fortran_type_ptr();
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    void* lidort_linout_lcl = lidort_linout_->fortran_type_ptr();

    lps_masters_write_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_Lps_Masters &obj) {
    output_stream << "Lidort_Lps_Masters:" << std::endl
      << "         thread: " << obj.thread()  << std::endl
      << "   lidort_fixin: " << obj.lidort_fixin()  << std::endl
      << "   lidort_modin: " << obj.lidort_modin()  << std::endl
      << "     lidort_sup: " << obj.lidort_sup()  << std::endl
      << "     lidort_out: " << obj.lidort_out()  << std::endl
      << "lidort_linfixin: " << obj.lidort_linfixin()  << std::endl
      << "lidort_linmodin: " << obj.lidort_linmodin()  << std::endl
      << "  lidort_linsup: " << obj.lidort_linsup()  << std::endl
      << "  lidort_linout: " << obj.lidort_linout()  << std::endl;
    return output_stream;
  }

private:
  int thread_;
  boost::shared_ptr<Lidort_Fixed_Inputs> lidort_fixin_;
  boost::shared_ptr<Lidort_Modified_Inputs> lidort_modin_;
  boost::shared_ptr<Lidort_Sup_Inout> lidort_sup_;
  boost::shared_ptr<Lidort_Outputs> lidort_out_;
  boost::shared_ptr<Lidort_Fixed_Lininputs> lidort_linfixin_;
  boost::shared_ptr<Lidort_Modified_Lininputs> lidort_linmodin_;
  boost::shared_ptr<Lidort_Linsup_Inout> lidort_linsup_;
  boost::shared_ptr<Lidort_Linoutputs> lidort_linout_;
};

//-----------------------------------------------------------------------
// Links to module: "lidort_inputs" in file: "lidort_inputs.F90"
//-----------------------------------------------------------------------

extern "C" {
  void inputs_read_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_inputstatus_in);
  void inputs_write_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_inputstatus_in);
  void inputs_input_master_wrap(const int* filnam_in_len, const char* filnam_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_inputstatus_in);
}

class Lidort_Inputs {

public:
  Lidort_Inputs() 
  { 
    
    // Initialize type pointers
    lidort_fixin_.reset( new Lidort_Fixed_Inputs() );
    lidort_modin_.reset( new Lidort_Modified_Inputs() );
    lidort_inputstatus_.reset( new Lidort_Input_Exception_Handling() );
    
  }

  Lidort_Fixed_Inputs& lidort_fixin() {
    return *lidort_fixin_;
  }
  const Lidort_Fixed_Inputs& lidort_fixin() const {
    return *lidort_fixin_;
  }
  boost::shared_ptr<Lidort_Fixed_Inputs>& lidort_fixin_ptr() {
    return lidort_fixin_;
  }
  

  Lidort_Modified_Inputs& lidort_modin() {
    return *lidort_modin_;
  }
  const Lidort_Modified_Inputs& lidort_modin() const {
    return *lidort_modin_;
  }
  boost::shared_ptr<Lidort_Modified_Inputs>& lidort_modin_ptr() {
    return lidort_modin_;
  }
  

  Lidort_Input_Exception_Handling& lidort_inputstatus() {
    return *lidort_inputstatus_;
  }
  const Lidort_Input_Exception_Handling& lidort_inputstatus() const {
    return *lidort_inputstatus_;
  }
  boost::shared_ptr<Lidort_Input_Exception_Handling>& lidort_inputstatus_ptr() {
    return lidort_inputstatus_;
  }
  

  
  void read_config(const std::string& filnam_in) {
    const char* filnam_lcl = filnam_in.c_str();
    int filnam_in_len = (int) filnam_in.size();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_inputstatus_lcl = lidort_inputstatus_->fortran_type_ptr();
    
    inputs_input_master_wrap(&filnam_in_len, filnam_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_inputstatus_lcl);
    

  }
 
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_inputstatus_lcl = lidort_inputstatus_->fortran_type_ptr();

    inputs_read_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_inputstatus_lcl);
    
  }
   
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_inputstatus_lcl = lidort_inputstatus_->fortran_type_ptr();

    inputs_write_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_inputstatus_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_Inputs &obj) {
    output_stream << "Lidort_Inputs:" << std::endl
      << "      lidort_fixin: " << obj.lidort_fixin()  << std::endl
      << "      lidort_modin: " << obj.lidort_modin()  << std::endl
      << "lidort_inputstatus: " << obj.lidort_inputstatus()  << std::endl;
    return output_stream;
  }

private:
  boost::shared_ptr<Lidort_Fixed_Inputs> lidort_fixin_;
  boost::shared_ptr<Lidort_Modified_Inputs> lidort_modin_;
  boost::shared_ptr<Lidort_Input_Exception_Handling> lidort_inputstatus_;
};

//-----------------------------------------------------------------------
// Links to module: "lidort_masters" in file: "lidort_masters.F90"
//-----------------------------------------------------------------------

extern "C" {
  void masters_read_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in);
  void masters_write_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in);
  void masters_master_wrap(const int* thread_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in);
}

class Lidort_Masters {

public:
  Lidort_Masters(const int& thread_in) : thread_(thread_in) 
  { 
    
    // Initialize type pointers
    lidort_fixin_.reset( new Lidort_Fixed_Inputs() );
    lidort_modin_.reset( new Lidort_Modified_Inputs() );
    lidort_sup_.reset( new Lidort_Sup_Inout() );
    lidort_out_.reset( new Lidort_Outputs() );
    
  }

  const int& thread() const {
    return thread_;
  }
  

  Lidort_Fixed_Inputs& lidort_fixin() {
    return *lidort_fixin_;
  }
  const Lidort_Fixed_Inputs& lidort_fixin() const {
    return *lidort_fixin_;
  }
  boost::shared_ptr<Lidort_Fixed_Inputs>& lidort_fixin_ptr() {
    return lidort_fixin_;
  }
  void lidort_fixin(Lidort_Fixed_Inputs& lidort_fixin_in) {
    void* src_ptr = lidort_fixin_in.fortran_type_ptr();
    void* dst_ptr = lidort_fixin_->fortran_type_ptr();
    lidort_fixed_inputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Modified_Inputs& lidort_modin() {
    return *lidort_modin_;
  }
  const Lidort_Modified_Inputs& lidort_modin() const {
    return *lidort_modin_;
  }
  boost::shared_ptr<Lidort_Modified_Inputs>& lidort_modin_ptr() {
    return lidort_modin_;
  }
  void lidort_modin(Lidort_Modified_Inputs& lidort_modin_in) {
    void* src_ptr = lidort_modin_in.fortran_type_ptr();
    void* dst_ptr = lidort_modin_->fortran_type_ptr();
    lidort_modified_inputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Sup_Inout& lidort_sup() {
    return *lidort_sup_;
  }
  const Lidort_Sup_Inout& lidort_sup() const {
    return *lidort_sup_;
  }
  boost::shared_ptr<Lidort_Sup_Inout>& lidort_sup_ptr() {
    return lidort_sup_;
  }
  void lidort_sup(Lidort_Sup_Inout& lidort_sup_in) {
    void* src_ptr = lidort_sup_in.fortran_type_ptr();
    void* dst_ptr = lidort_sup_->fortran_type_ptr();
    lidort_sup_inout_c_copy(&src_ptr, &dst_ptr);
  }
  

  Lidort_Outputs& lidort_out() {
    return *lidort_out_;
  }
  const Lidort_Outputs& lidort_out() const {
    return *lidort_out_;
  }
  boost::shared_ptr<Lidort_Outputs>& lidort_out_ptr() {
    return lidort_out_;
  }
  void lidort_out(Lidort_Outputs& lidort_out_in) {
    void* src_ptr = lidort_out_in.fortran_type_ptr();
    void* dst_ptr = lidort_out_->fortran_type_ptr();
    lidort_outputs_c_copy(&src_ptr, &dst_ptr);
  }
  

  
  void run() {
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();
    
    masters_master_wrap(&thread_, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl);
    
    Lidort_Pars lid_pars = Lidort_Pars::instance();
    if( lidort_out().status().ts_status_inputcheck() != lid_pars.lidort_success ||
        lidort_out().status().ts_status_calculation() != lid_pars.lidort_success ) {
       std::stringstream err_msg;
       err_msg << "LIDORT Error at " << __FILE__ << ":" << __LINE__ << std::endl;
       // Output the full details of the error message to stderr since the exception
       // class may truncate the message
       std::cerr << err_msg.str();
       std::cerr << lidort_out().status();
       throw Exception(err_msg.str());
    }

  }
 
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();

    masters_read_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl);
    
  }
   
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();

    masters_write_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_Masters &obj) {
    output_stream << "Lidort_Masters:" << std::endl
      << "      thread: " << obj.thread()  << std::endl
      << "lidort_fixin: " << obj.lidort_fixin()  << std::endl
      << "lidort_modin: " << obj.lidort_modin()  << std::endl
      << "  lidort_sup: " << obj.lidort_sup()  << std::endl
      << "  lidort_out: " << obj.lidort_out()  << std::endl;
    return output_stream;
  }

private:
  int thread_;
  boost::shared_ptr<Lidort_Fixed_Inputs> lidort_fixin_;
  boost::shared_ptr<Lidort_Modified_Inputs> lidort_modin_;
  boost::shared_ptr<Lidort_Sup_Inout> lidort_sup_;
  boost::shared_ptr<Lidort_Outputs> lidort_out_;
};

//-----------------------------------------------------------------------
// Links to module: "lidort_sup_accessories" in file: "lidort_sup_accessories.F90"
//-----------------------------------------------------------------------

extern "C" {
  void sup_accessories_read_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_brdfcheck_status_in);
  void sup_accessories_write_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_brdfcheck_status_in);
  void sup_accessories_brdf_input_checker_wrap(void** brdf_sup_in_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_brdfcheck_status_in);
}

class Lidort_Sup_Accessories {

public:
  Lidort_Sup_Accessories(boost::shared_ptr<Brdf_Sup_Inputs>& brdf_sup_in_in, boost::shared_ptr<Lidort_Fixed_Inputs>& lidort_fixin_in, boost::shared_ptr<Lidort_Modified_Inputs>& lidort_modin_in) : brdf_sup_in_(brdf_sup_in_in), lidort_fixin_(lidort_fixin_in), lidort_modin_(lidort_modin_in) 
  { 
    
    // Initialize type pointers
    lidort_brdfcheck_status_.reset( new Lidort_Exception_Handling() );
    
  }

  Brdf_Sup_Inputs& brdf_sup_in() {
    return *brdf_sup_in_;
  }
  const Brdf_Sup_Inputs& brdf_sup_in() const {
    return *brdf_sup_in_;
  }
  boost::shared_ptr<Brdf_Sup_Inputs>& brdf_sup_in_ptr() {
    return brdf_sup_in_;
  }
  

  Lidort_Fixed_Inputs& lidort_fixin() {
    return *lidort_fixin_;
  }
  const Lidort_Fixed_Inputs& lidort_fixin() const {
    return *lidort_fixin_;
  }
  boost::shared_ptr<Lidort_Fixed_Inputs>& lidort_fixin_ptr() {
    return lidort_fixin_;
  }
  

  Lidort_Modified_Inputs& lidort_modin() {
    return *lidort_modin_;
  }
  const Lidort_Modified_Inputs& lidort_modin() const {
    return *lidort_modin_;
  }
  boost::shared_ptr<Lidort_Modified_Inputs>& lidort_modin_ptr() {
    return lidort_modin_;
  }
  

  Lidort_Exception_Handling& lidort_brdfcheck_status() {
    return *lidort_brdfcheck_status_;
  }
  const Lidort_Exception_Handling& lidort_brdfcheck_status() const {
    return *lidort_brdfcheck_status_;
  }
  boost::shared_ptr<Lidort_Exception_Handling>& lidort_brdfcheck_status_ptr() {
    return lidort_brdfcheck_status_;
  }
  

  
  void brdf_input_checker() {
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_brdfcheck_status_lcl = lidort_brdfcheck_status_->fortran_type_ptr();
    
    sup_accessories_brdf_input_checker_wrap(&brdf_sup_in_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_brdfcheck_status_lcl);
    

  }
 
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_brdfcheck_status_lcl = lidort_brdfcheck_status_->fortran_type_ptr();

    sup_accessories_read_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_brdfcheck_status_lcl);
    
  }
   
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_brdfcheck_status_lcl = lidort_brdfcheck_status_->fortran_type_ptr();

    sup_accessories_write_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_brdfcheck_status_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_Sup_Accessories &obj) {
    output_stream << "Lidort_Sup_Accessories:" << std::endl
      << "            brdf_sup_in: " << obj.brdf_sup_in()  << std::endl
      << "           lidort_fixin: " << obj.lidort_fixin()  << std::endl
      << "           lidort_modin: " << obj.lidort_modin()  << std::endl
      << "lidort_brdfcheck_status: " << obj.lidort_brdfcheck_status()  << std::endl;
    return output_stream;
  }

private:
  boost::shared_ptr<Brdf_Sup_Inputs> brdf_sup_in_;
  boost::shared_ptr<Lidort_Fixed_Inputs> lidort_fixin_;
  boost::shared_ptr<Lidort_Modified_Inputs> lidort_modin_;
  boost::shared_ptr<Lidort_Exception_Handling> lidort_brdfcheck_status_;
};



}
#endif