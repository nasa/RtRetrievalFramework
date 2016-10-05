#ifndef GLOBAL_FIXTURE_H
#define GLOBAL_FIXTURE_H
#include <string>
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This is a global fixture that is available to all unit tests.
*******************************************************************/
class GlobalFixture {
public:
  GlobalFixture();
  virtual ~GlobalFixture();
  std::string input_dir() const;
  std::string test_data_dir() const;
  std::string absco_data_dir() const;
  std::string merra_data_dir() const;
  std::string absco_4d_dir() const;
  void turn_on_logger() const;
  // Add file to be removed at end of unit test. Ok if file doesn't
  // actually exist (i.e., no special handling needed for error
  // conditions). 
  void add_file_to_cleanup(const std::string& Fname)
  { cleanup_list.push_back(Fname); }
private:
  void set_default_value();
  std::vector<std::string> cleanup_list;
};
}
#endif
