#ifndef LIDORT_INTERFACE_MASTERS_H
#define LIDORT_INTERFACE_MASTERS_H

#include <iostream>
#include <blitz/array.h>

#include "fp_exception.h"
#include "lidort_interface_types.h"

/* This file was auto-generated */

namespace FullPhysics {

//-----------------------------------------------------------------------
// Links to module: "brdf_lin_sup_masters_m" in file: "brdf_lin_sup_masters.F90"
//-----------------------------------------------------------------------

extern "C" {
  void brdf_lin_sup_masters_m_read_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** brdf_linsup_in_in, void** brdf_sup_inputstatus_in, void** brdf_sup_out_in, void** brdf_linsup_out_in, void** brdf_sup_outputstatus_in);
  void brdf_lin_sup_masters_m_write_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** brdf_linsup_in_in, void** brdf_sup_inputstatus_in, void** brdf_sup_out_in, void** brdf_linsup_out_in, void** brdf_sup_outputstatus_in);
  void brdf_lin_sup_masters_m_brdf_lin_inputmaster_wrap(const int* filnam_in_len, const char* filnam_in, void** brdf_sup_in_in, void** brdf_linsup_in_in, void** brdf_sup_inputstatus_in);
  void brdf_lin_sup_masters_m_brdf_lin_mainmaster_wrap(const bool* do_debug_restoration_in, const int* nmoments_input_in, void** brdf_sup_in_in, void** brdf_linsup_in_in, void** brdf_sup_out_in, void** brdf_linsup_out_in, void** brdf_sup_outputstatus_in);
}

class Brdf_Lin_Sup_Masters : public virtual GenericObject {

public:
  Brdf_Lin_Sup_Masters() 
  { 
    
    // Initialize type pointers
    brdf_sup_in_.reset( new Brdf_Sup_Inputs() );
    brdf_linsup_in_.reset( new Brdf_Linsup_Inputs() );
    brdf_sup_inputstatus_.reset( new Brdf_Input_Exception_Handling() );
    brdf_sup_out_.reset( new Brdf_Sup_Outputs() );
    brdf_linsup_out_.reset( new Brdf_Linsup_Outputs() );
    brdf_sup_outputstatus_.reset( new Brdf_Output_Exception_Handling() );
    
  }

  virtual ~Brdf_Lin_Sup_Masters() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
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

  

  Brdf_Output_Exception_Handling& brdf_sup_outputstatus() {
    return *brdf_sup_outputstatus_;
  }

  const Brdf_Output_Exception_Handling& brdf_sup_outputstatus() const {
    return *brdf_sup_outputstatus_;
  }

  boost::shared_ptr<Brdf_Output_Exception_Handling>& brdf_sup_outputstatus_ptr() {
    return brdf_sup_outputstatus_;
  }

  

  
  void read_config(const std::string& filnam_in) {
    const char* filnam_lcl = filnam_in.c_str();
    int filnam_in_len = (int) filnam_in.size();
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_linsup_in_lcl = brdf_linsup_in_->fortran_type_ptr();
    void* brdf_sup_inputstatus_lcl = brdf_sup_inputstatus_->fortran_type_ptr();
    
    brdf_lin_sup_masters_m_brdf_lin_inputmaster_wrap(&filnam_in_len, filnam_lcl, &brdf_sup_in_lcl, &brdf_linsup_in_lcl, &brdf_sup_inputstatus_lcl);
    

  }
void run(const bool& do_debug_restoration_in, const int& nmoments_input_in) {
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_linsup_in_lcl = brdf_linsup_in_->fortran_type_ptr();
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();
    void* brdf_linsup_out_lcl = brdf_linsup_out_->fortran_type_ptr();
    void* brdf_sup_outputstatus_lcl = brdf_sup_outputstatus_->fortran_type_ptr();
    
    brdf_lin_sup_masters_m_brdf_lin_mainmaster_wrap(&do_debug_restoration_in, &nmoments_input_in, &brdf_sup_in_lcl, &brdf_linsup_in_lcl, &brdf_sup_out_lcl, &brdf_linsup_out_lcl, &brdf_sup_outputstatus_lcl);
    

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
    void* brdf_sup_outputstatus_lcl = brdf_sup_outputstatus_->fortran_type_ptr();

    brdf_lin_sup_masters_m_read_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &brdf_linsup_in_lcl, &brdf_sup_inputstatus_lcl, &brdf_sup_out_lcl, &brdf_linsup_out_lcl, &brdf_sup_outputstatus_lcl);
    
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
    void* brdf_sup_outputstatus_lcl = brdf_sup_outputstatus_->fortran_type_ptr();

    brdf_lin_sup_masters_m_write_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &brdf_linsup_in_lcl, &brdf_sup_inputstatus_lcl, &brdf_sup_out_lcl, &brdf_linsup_out_lcl, &brdf_sup_outputstatus_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Brdf_Lin_Sup_Masters &obj) {
    output_stream << "Brdf_Lin_Sup_Masters:" << std::endl
      << "          brdf_sup_in: " << obj.brdf_sup_in()  << std::endl
      << "       brdf_linsup_in: " << obj.brdf_linsup_in()  << std::endl
      << " brdf_sup_inputstatus: " << obj.brdf_sup_inputstatus()  << std::endl
      << "         brdf_sup_out: " << obj.brdf_sup_out()  << std::endl
      << "      brdf_linsup_out: " << obj.brdf_linsup_out()  << std::endl
      << "brdf_sup_outputstatus: " << obj.brdf_sup_outputstatus()  << std::endl;
    return output_stream;

  }

private:
  boost::shared_ptr<Brdf_Sup_Inputs> brdf_sup_in_;
  boost::shared_ptr<Brdf_Linsup_Inputs> brdf_linsup_in_;
  boost::shared_ptr<Brdf_Input_Exception_Handling> brdf_sup_inputstatus_;
  boost::shared_ptr<Brdf_Sup_Outputs> brdf_sup_out_;
  boost::shared_ptr<Brdf_Linsup_Outputs> brdf_linsup_out_;
  boost::shared_ptr<Brdf_Output_Exception_Handling> brdf_sup_outputstatus_;
};

//-----------------------------------------------------------------------
// Links to module: "brdf_sup_masters_m" in file: "brdf_sup_masters.F90"
//-----------------------------------------------------------------------

extern "C" {
  void brdf_sup_masters_m_read_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** brdf_sup_inputstatus_in, void** brdf_sup_out_in, void** brdf_sup_outputstatus_in);
  void brdf_sup_masters_m_write_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** brdf_sup_inputstatus_in, void** brdf_sup_out_in, void** brdf_sup_outputstatus_in);
  void brdf_sup_masters_m_brdf_inputmaster_wrap(const int* filnam_in_len, const char* filnam_in, void** brdf_sup_in_in, void** brdf_sup_inputstatus_in);
  void brdf_sup_masters_m_brdf_mainmaster_wrap(const bool* do_debug_restoration_in, const int* nmoments_input_in, void** brdf_sup_in_in, void** brdf_sup_out_in, void** brdf_sup_outputstatus_in);
}

class Brdf_Sup_Masters : public virtual GenericObject {

public:
  Brdf_Sup_Masters() 
  { 
    
    // Initialize type pointers
    brdf_sup_in_.reset( new Brdf_Sup_Inputs() );
    brdf_sup_inputstatus_.reset( new Brdf_Input_Exception_Handling() );
    brdf_sup_out_.reset( new Brdf_Sup_Outputs() );
    brdf_sup_outputstatus_.reset( new Brdf_Output_Exception_Handling() );
    
  }

  virtual ~Brdf_Sup_Masters() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
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

  

  Brdf_Sup_Outputs& brdf_sup_out() {
    return *brdf_sup_out_;
  }

  const Brdf_Sup_Outputs& brdf_sup_out() const {
    return *brdf_sup_out_;
  }

  boost::shared_ptr<Brdf_Sup_Outputs>& brdf_sup_out_ptr() {
    return brdf_sup_out_;
  }

  

  Brdf_Output_Exception_Handling& brdf_sup_outputstatus() {
    return *brdf_sup_outputstatus_;
  }

  const Brdf_Output_Exception_Handling& brdf_sup_outputstatus() const {
    return *brdf_sup_outputstatus_;
  }

  boost::shared_ptr<Brdf_Output_Exception_Handling>& brdf_sup_outputstatus_ptr() {
    return brdf_sup_outputstatus_;
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
    void* brdf_sup_outputstatus_lcl = brdf_sup_outputstatus_->fortran_type_ptr();
    
    brdf_sup_masters_m_brdf_mainmaster_wrap(&do_debug_restoration_in, &nmoments_input_in, &brdf_sup_in_lcl, &brdf_sup_out_lcl, &brdf_sup_outputstatus_lcl);
    

  }
 
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_sup_inputstatus_lcl = brdf_sup_inputstatus_->fortran_type_ptr();
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();
    void* brdf_sup_outputstatus_lcl = brdf_sup_outputstatus_->fortran_type_ptr();

    brdf_sup_masters_m_read_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &brdf_sup_inputstatus_lcl, &brdf_sup_out_lcl, &brdf_sup_outputstatus_lcl);
    
  }
   
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* brdf_sup_inputstatus_lcl = brdf_sup_inputstatus_->fortran_type_ptr();
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();
    void* brdf_sup_outputstatus_lcl = brdf_sup_outputstatus_->fortran_type_ptr();

    brdf_sup_masters_m_write_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &brdf_sup_inputstatus_lcl, &brdf_sup_out_lcl, &brdf_sup_outputstatus_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Brdf_Sup_Masters &obj) {
    output_stream << "Brdf_Sup_Masters:" << std::endl
      << "          brdf_sup_in: " << obj.brdf_sup_in()  << std::endl
      << " brdf_sup_inputstatus: " << obj.brdf_sup_inputstatus()  << std::endl
      << "         brdf_sup_out: " << obj.brdf_sup_out()  << std::endl
      << "brdf_sup_outputstatus: " << obj.brdf_sup_outputstatus()  << std::endl;
    return output_stream;

  }

private:
  boost::shared_ptr<Brdf_Sup_Inputs> brdf_sup_in_;
  boost::shared_ptr<Brdf_Input_Exception_Handling> brdf_sup_inputstatus_;
  boost::shared_ptr<Brdf_Sup_Outputs> brdf_sup_out_;
  boost::shared_ptr<Brdf_Output_Exception_Handling> brdf_sup_outputstatus_;
};

//-----------------------------------------------------------------------
// Links to module: "lidort_inputs_m" in file: "lidort_inputs.F90"
//-----------------------------------------------------------------------

extern "C" {
  void inputs_m_read_wrap(const char* filename, const int* filename_len, void** lidort_sup_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_inputstatus_in);
  void inputs_m_write_wrap(const char* filename, const int* filename_len, void** lidort_sup_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_inputstatus_in);
  void inputs_m_brdf_sup_init_wrap(void** lidort_sup_in);
  void inputs_m_input_master_wrap(const int* filnam_in_len, const char* filnam_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_inputstatus_in);
  void inputs_m_sleave_sup_init_wrap(void** lidort_sup_in);
  void inputs_m_ss_sup_init_wrap(void** lidort_sup_in);
  void inputs_m_sup_init_wrap(void** lidort_sup_in);
}

class Lidort_Inputs : public virtual GenericObject {

public:
  Lidort_Inputs() 
  { 
    
    // Initialize type pointers
    lidort_sup_.reset( new Lidort_Sup_Inout() );
    lidort_fixin_.reset( new Lidort_Fixed_Inputs() );
    lidort_modin_.reset( new Lidort_Modified_Inputs() );
    lidort_inputstatus_.reset( new Lidort_Input_Exception_Handling() );
    
  }

  virtual ~Lidort_Inputs() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
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

  

  
  void brdf_sup_init() {
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    
    inputs_m_brdf_sup_init_wrap(&lidort_sup_lcl);
    

  }
void read_config(const std::string& filnam_in) {
    const char* filnam_lcl = filnam_in.c_str();
    int filnam_in_len = (int) filnam_in.size();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_inputstatus_lcl = lidort_inputstatus_->fortran_type_ptr();
    
    inputs_m_input_master_wrap(&filnam_in_len, filnam_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_inputstatus_lcl);
    

  }
void sleave_sup_init() {
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    
    inputs_m_sleave_sup_init_wrap(&lidort_sup_lcl);
    

  }
void ss_sup_init() {
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    
    inputs_m_ss_sup_init_wrap(&lidort_sup_lcl);
    

  }
void sup_init() {
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    
    inputs_m_sup_init_wrap(&lidort_sup_lcl);
    

  }
 
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_inputstatus_lcl = lidort_inputstatus_->fortran_type_ptr();

    inputs_m_read_wrap(filename_lcl, &filename_in_len, &lidort_sup_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_inputstatus_lcl);
    
  }
   
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_inputstatus_lcl = lidort_inputstatus_->fortran_type_ptr();

    inputs_m_write_wrap(filename_lcl, &filename_in_len, &lidort_sup_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_inputstatus_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_Inputs &obj) {
    output_stream << "Lidort_Inputs:" << std::endl
      << "        lidort_sup: " << obj.lidort_sup()  << std::endl
      << "      lidort_fixin: " << obj.lidort_fixin()  << std::endl
      << "      lidort_modin: " << obj.lidort_modin()  << std::endl
      << "lidort_inputstatus: " << obj.lidort_inputstatus()  << std::endl;
    return output_stream;

  }

private:
  boost::shared_ptr<Lidort_Sup_Inout> lidort_sup_;
  boost::shared_ptr<Lidort_Fixed_Inputs> lidort_fixin_;
  boost::shared_ptr<Lidort_Modified_Inputs> lidort_modin_;
  boost::shared_ptr<Lidort_Input_Exception_Handling> lidort_inputstatus_;
};

//-----------------------------------------------------------------------
// Links to module: "lidort_masters_m" in file: "lidort_masters.F90"
//-----------------------------------------------------------------------

extern "C" {
  void masters_m_read_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in);
  void masters_m_write_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in);
  void masters_m_master_wrap(const bool* do_debug_input_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in);
}

class Lidort_Masters : public virtual GenericObject {

public:
  Lidort_Masters() 
  { 
    
    // Initialize type pointers
    lidort_fixin_.reset( new Lidort_Fixed_Inputs() );
    lidort_modin_.reset( new Lidort_Modified_Inputs() );
    lidort_sup_.reset( new Lidort_Sup_Inout() );
    lidort_out_.reset( new Lidort_Outputs() );
    
  }

  virtual ~Lidort_Masters() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
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

  

  
  void run(const bool& do_debug_input_in) {
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();
    
    masters_m_master_wrap(&do_debug_input_in, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl);
    
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

    masters_m_read_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl);
    
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

    masters_m_write_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_Masters &obj) {
    output_stream << "Lidort_Masters:" << std::endl
      << "lidort_fixin: " << obj.lidort_fixin()  << std::endl
      << "lidort_modin: " << obj.lidort_modin()  << std::endl
      << "  lidort_sup: " << obj.lidort_sup()  << std::endl
      << "  lidort_out: " << obj.lidort_out()  << std::endl;
    return output_stream;

  }

private:
  boost::shared_ptr<Lidort_Fixed_Inputs> lidort_fixin_;
  boost::shared_ptr<Lidort_Modified_Inputs> lidort_modin_;
  boost::shared_ptr<Lidort_Sup_Inout> lidort_sup_;
  boost::shared_ptr<Lidort_Outputs> lidort_out_;
};

//-----------------------------------------------------------------------
// Links to module: "lidort_l_inputs_m" in file: "lidort_l_inputs.F90"
//-----------------------------------------------------------------------

extern "C" {
  void l_inputs_m_read_wrap(const char* filename, const int* filename_len, void** lidort_linsup_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_inputstatus_in);
  void l_inputs_m_write_wrap(const char* filename, const int* filename_len, void** lidort_linsup_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_inputstatus_in);
  void l_inputs_m_brdf_linsup_init_wrap(void** lidort_linsup_in);
  void l_inputs_m_l_input_master_wrap(const int* filnam_in_len, const char* filnam_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_inputstatus_in);
  void l_inputs_m_linsup_init_wrap(void** lidort_linsup_in);
  void l_inputs_m_sleave_linsup_init_wrap(void** lidort_linsup_in);
  void l_inputs_m_ss_linsup_init_wrap(void** lidort_linsup_in);
}

class Lidort_L_Inputs : public virtual GenericObject {

public:
  Lidort_L_Inputs() 
  { 
    
    // Initialize type pointers
    lidort_linsup_.reset( new Lidort_Linsup_Inout() );
    lidort_fixin_.reset( new Lidort_Fixed_Inputs() );
    lidort_modin_.reset( new Lidort_Modified_Inputs() );
    lidort_linfixin_.reset( new Lidort_Fixed_Lininputs() );
    lidort_linmodin_.reset( new Lidort_Modified_Lininputs() );
    lidort_inputstatus_.reset( new Lidort_Input_Exception_Handling() );
    
  }

  virtual ~Lidort_L_Inputs() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
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

  

  Lidort_Fixed_Lininputs& lidort_linfixin() {
    return *lidort_linfixin_;
  }

  const Lidort_Fixed_Lininputs& lidort_linfixin() const {
    return *lidort_linfixin_;
  }

  boost::shared_ptr<Lidort_Fixed_Lininputs>& lidort_linfixin_ptr() {
    return lidort_linfixin_;
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

  

  Lidort_Input_Exception_Handling& lidort_inputstatus() {
    return *lidort_inputstatus_;
  }

  const Lidort_Input_Exception_Handling& lidort_inputstatus() const {
    return *lidort_inputstatus_;
  }

  boost::shared_ptr<Lidort_Input_Exception_Handling>& lidort_inputstatus_ptr() {
    return lidort_inputstatus_;
  }

  

  
  void brdf_linsup_init() {
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    
    l_inputs_m_brdf_linsup_init_wrap(&lidort_linsup_lcl);
    

  }
void read_config(const std::string& filnam_in) {
    const char* filnam_lcl = filnam_in.c_str();
    int filnam_in_len = (int) filnam_in.size();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_linfixin_lcl = lidort_linfixin_->fortran_type_ptr();
    void* lidort_linmodin_lcl = lidort_linmodin_->fortran_type_ptr();
    void* lidort_inputstatus_lcl = lidort_inputstatus_->fortran_type_ptr();
    
    l_inputs_m_l_input_master_wrap(&filnam_in_len, filnam_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_inputstatus_lcl);
    

  }
void linsup_init() {
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    
    l_inputs_m_linsup_init_wrap(&lidort_linsup_lcl);
    

  }
void sleave_linsup_init() {
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    
    l_inputs_m_sleave_linsup_init_wrap(&lidort_linsup_lcl);
    

  }
void ss_linsup_init() {
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    
    l_inputs_m_ss_linsup_init_wrap(&lidort_linsup_lcl);
    

  }
 
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void read_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_linfixin_lcl = lidort_linfixin_->fortran_type_ptr();
    void* lidort_linmodin_lcl = lidort_linmodin_->fortran_type_ptr();
    void* lidort_inputstatus_lcl = lidort_inputstatus_->fortran_type_ptr();

    l_inputs_m_read_wrap(filename_lcl, &filename_in_len, &lidort_linsup_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_inputstatus_lcl);
    
  }
   
  // This routine is meant only for testing purposes and interoperability 
  // with offline tests
  void write_fortran_file(const std::string& filename_in) {
    const char* filename_lcl = filename_in.c_str();
    int filename_in_len = (int) filename_in.size();
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_linfixin_lcl = lidort_linfixin_->fortran_type_ptr();
    void* lidort_linmodin_lcl = lidort_linmodin_->fortran_type_ptr();
    void* lidort_inputstatus_lcl = lidort_inputstatus_->fortran_type_ptr();

    l_inputs_m_write_wrap(filename_lcl, &filename_in_len, &lidort_linsup_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_inputstatus_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_L_Inputs &obj) {
    output_stream << "Lidort_L_Inputs:" << std::endl
      << "     lidort_linsup: " << obj.lidort_linsup()  << std::endl
      << "      lidort_fixin: " << obj.lidort_fixin()  << std::endl
      << "      lidort_modin: " << obj.lidort_modin()  << std::endl
      << "   lidort_linfixin: " << obj.lidort_linfixin()  << std::endl
      << "   lidort_linmodin: " << obj.lidort_linmodin()  << std::endl
      << "lidort_inputstatus: " << obj.lidort_inputstatus()  << std::endl;
    return output_stream;

  }

private:
  boost::shared_ptr<Lidort_Linsup_Inout> lidort_linsup_;
  boost::shared_ptr<Lidort_Fixed_Inputs> lidort_fixin_;
  boost::shared_ptr<Lidort_Modified_Inputs> lidort_modin_;
  boost::shared_ptr<Lidort_Fixed_Lininputs> lidort_linfixin_;
  boost::shared_ptr<Lidort_Modified_Lininputs> lidort_linmodin_;
  boost::shared_ptr<Lidort_Input_Exception_Handling> lidort_inputstatus_;
};

//-----------------------------------------------------------------------
// Links to module: "lidort_lcs_masters_m" in file: "lidort_lcs_masters.F90"
//-----------------------------------------------------------------------

extern "C" {
  void lcs_masters_m_read_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
  void lcs_masters_m_write_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
  void lcs_masters_m_lcs_master_wrap(const bool* do_debug_input_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
}

class Lidort_Lcs_Masters : public virtual GenericObject {

public:
  Lidort_Lcs_Masters() 
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

  virtual ~Lidort_Lcs_Masters() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
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

  

  
  void run(const bool& do_debug_input_in) {
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();
    void* lidort_linfixin_lcl = lidort_linfixin_->fortran_type_ptr();
    void* lidort_linmodin_lcl = lidort_linmodin_->fortran_type_ptr();
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    void* lidort_linout_lcl = lidort_linout_->fortran_type_ptr();
    
    lcs_masters_m_lcs_master_wrap(&do_debug_input_in, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
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

    lcs_masters_m_read_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
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

    lcs_masters_m_write_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_Lcs_Masters &obj) {
    output_stream << "Lidort_Lcs_Masters:" << std::endl
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
// Links to module: "lidort_lps_masters_m" in file: "lidort_lps_masters.F90"
//-----------------------------------------------------------------------

extern "C" {
  void lps_masters_m_read_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
  void lps_masters_m_write_wrap(const char* filename, const int* filename_len, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
  void lps_masters_m_lps_master_wrap(const bool* do_debug_input_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in, void** lidort_out_in, void** lidort_linfixin_in, void** lidort_linmodin_in, void** lidort_linsup_in, void** lidort_linout_in);
}

class Lidort_Lps_Masters : public virtual GenericObject {

public:
  Lidort_Lps_Masters() 
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

  virtual ~Lidort_Lps_Masters() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
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

  

  
  void run(const bool& do_debug_input_in) {
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    void* lidort_out_lcl = lidort_out_->fortran_type_ptr();
    void* lidort_linfixin_lcl = lidort_linfixin_->fortran_type_ptr();
    void* lidort_linmodin_lcl = lidort_linmodin_->fortran_type_ptr();
    void* lidort_linsup_lcl = lidort_linsup_->fortran_type_ptr();
    void* lidort_linout_lcl = lidort_linout_->fortran_type_ptr();
    
    lps_masters_m_lps_master_wrap(&do_debug_input_in, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
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

    lps_masters_m_read_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
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

    lps_masters_m_write_wrap(filename_lcl, &filename_in_len, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl, &lidort_out_lcl, &lidort_linfixin_lcl, &lidort_linmodin_lcl, &lidort_linsup_lcl, &lidort_linout_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_Lps_Masters &obj) {
    output_stream << "Lidort_Lps_Masters:" << std::endl
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
// Links to module: "lidort_brdf_sup_accessories_m" in file: "lidort_brdf_sup_accessories.F90"
//-----------------------------------------------------------------------

extern "C" {
  void brdf_sup_accessories_m_read_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_brdfcheck_status_in, void** brdf_sup_out_in, void** lidort_sup_in);
  void brdf_sup_accessories_m_write_wrap(const char* filename, const int* filename_len, void** brdf_sup_in_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_brdfcheck_status_in, void** brdf_sup_out_in, void** lidort_sup_in);
  void brdf_sup_accessories_m_brdf_input_check_wrap(void** brdf_sup_in_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_brdfcheck_status_in);
  void brdf_sup_accessories_m_brdf_input_check_error_wrap(const int* errorfile_in_len, const char* errorfile_in, void** lidort_brdfcheck_status_in);
  void brdf_sup_accessories_m_set_brdf_inputs_wrap(void** brdf_sup_out_in, void** lidort_fixin_in, void** lidort_modin_in, void** lidort_sup_in);
}

class Lidort_Brdf_Sup_Accessories : public virtual GenericObject {

public:
  Lidort_Brdf_Sup_Accessories(boost::shared_ptr<Brdf_Sup_Inputs>& brdf_sup_in_in, boost::shared_ptr<Lidort_Fixed_Inputs>& lidort_fixin_in, boost::shared_ptr<Lidort_Modified_Inputs>& lidort_modin_in) : brdf_sup_in_(brdf_sup_in_in), lidort_fixin_(lidort_fixin_in), lidort_modin_(lidort_modin_in) 
  { 
    
    // Initialize type pointers
    lidort_brdfcheck_status_.reset( new Lidort_Exception_Handling() );
    brdf_sup_out_.reset( new Brdf_Sup_Outputs() );
    lidort_sup_.reset( new Lidort_Sup_Inout() );
    
  }

  virtual ~Lidort_Brdf_Sup_Accessories() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
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

  

  Brdf_Sup_Outputs& brdf_sup_out() {
    return *brdf_sup_out_;
  }

  const Brdf_Sup_Outputs& brdf_sup_out() const {
    return *brdf_sup_out_;
  }

  boost::shared_ptr<Brdf_Sup_Outputs>& brdf_sup_out_ptr() {
    return brdf_sup_out_;
  }

  void brdf_sup_out(Brdf_Sup_Outputs& brdf_sup_out_in) {
    void* src_ptr = brdf_sup_out_in.fortran_type_ptr();
    void* dst_ptr = brdf_sup_out_->fortran_type_ptr();
    brdf_sup_outputs_c_copy(&src_ptr, &dst_ptr);
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

  

  
  void brdf_input_check() {
    void* brdf_sup_in_lcl = brdf_sup_in_->fortran_type_ptr();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_brdfcheck_status_lcl = lidort_brdfcheck_status_->fortran_type_ptr();
    
    brdf_sup_accessories_m_brdf_input_check_wrap(&brdf_sup_in_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_brdfcheck_status_lcl);
    

  }
void brdf_input_check_error(const std::string& errorfile_in) {
    const char* errorfile_lcl = errorfile_in.c_str();
    int errorfile_in_len = (int) errorfile_in.size();
    void* lidort_brdfcheck_status_lcl = lidort_brdfcheck_status_->fortran_type_ptr();
    
    brdf_sup_accessories_m_brdf_input_check_error_wrap(&errorfile_in_len, errorfile_lcl, &lidort_brdfcheck_status_lcl);
    

  }
void set_brdf_inputs() {
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();
    void* lidort_fixin_lcl = lidort_fixin_->fortran_type_ptr();
    void* lidort_modin_lcl = lidort_modin_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();
    
    brdf_sup_accessories_m_set_brdf_inputs_wrap(&brdf_sup_out_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_sup_lcl);
    

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
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();

    brdf_sup_accessories_m_read_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_brdfcheck_status_lcl, &brdf_sup_out_lcl, &lidort_sup_lcl);
    
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
    void* brdf_sup_out_lcl = brdf_sup_out_->fortran_type_ptr();
    void* lidort_sup_lcl = lidort_sup_->fortran_type_ptr();

    brdf_sup_accessories_m_write_wrap(filename_lcl, &filename_in_len, &brdf_sup_in_lcl, &lidort_fixin_lcl, &lidort_modin_lcl, &lidort_brdfcheck_status_lcl, &brdf_sup_out_lcl, &lidort_sup_lcl);
    
  }
  
  friend std::ostream& operator<<(std::ostream &output_stream, const Lidort_Brdf_Sup_Accessories &obj) {
    output_stream << "Lidort_Brdf_Sup_Accessories:" << std::endl
      << "            brdf_sup_in: " << obj.brdf_sup_in()  << std::endl
      << "           lidort_fixin: " << obj.lidort_fixin()  << std::endl
      << "           lidort_modin: " << obj.lidort_modin()  << std::endl
      << "lidort_brdfcheck_status: " << obj.lidort_brdfcheck_status()  << std::endl
      << "           brdf_sup_out: " << obj.brdf_sup_out()  << std::endl
      << "             lidort_sup: " << obj.lidort_sup()  << std::endl;
    return output_stream;

  }

private:
  boost::shared_ptr<Brdf_Sup_Inputs> brdf_sup_in_;
  boost::shared_ptr<Lidort_Fixed_Inputs> lidort_fixin_;
  boost::shared_ptr<Lidort_Modified_Inputs> lidort_modin_;
  boost::shared_ptr<Lidort_Exception_Handling> lidort_brdfcheck_status_;
  boost::shared_ptr<Brdf_Sup_Outputs> brdf_sup_out_;
  boost::shared_ptr<Lidort_Sup_Inout> lidort_sup_;
};



}
#endif