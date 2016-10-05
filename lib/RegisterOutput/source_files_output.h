#ifndef SOURCE_FILES_H
#define SOURCE_FILES_H

#include "register_output_base.h"
#include <vector>

namespace FullPhysics {

/****************************************************************//**
  Writes source filenames into the output file
*******************************************************************/
class SourceFilesOutput : public RegisterOutputBase {
public:
  SourceFilesOutput(const std::string& group_name, const std::vector<std::string>& dataset_names, const std::vector<std::string>& file_names) : 
    group_name_(group_name), dataset_names_(dataset_names), file_names_(file_names) {}
  virtual ~SourceFilesOutput() {}
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
private:
  std::string group_name_;
  std::vector<std::string> dataset_names_;
  std::vector<std::string> file_names_;
};
}
#endif
