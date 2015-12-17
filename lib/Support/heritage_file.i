// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "heritage_file.h"
%}

%base_import(generic_object)
%fp_shared_ptr(FullPhysics::HeritageFile);

namespace FullPhysics {
class HeritageFile : public GenericObject {
public:
  HeritageFile(const std::string& Fname);
  std::string print_to_string() const;
  void parse_file(const std::string& Fname);
  %python_attribute(data, const blitz::Array<double, 2>&);
  int column_index(const std::string& Col_name) const;
  blitz::Array<double, 1> data(const std::string& Col_name) const;
  bool has_value(const std::string& Keyword) const;
  %extend {
    int value_int(const std::string& Keyword) const {
       return $self->value<int>(Keyword);
    }
    double value_double(const std::string& Keyword) const {
       return $self->value<double>(Keyword);
    }
    std::string value_string(const std::string& Keyword) const {
       return $self->value<std::string>(Keyword);
    }
    bool value_bool(const std::string& Keyword) const {
       return $self->value<bool>(Keyword);
    }
    std::vector<std::string> value_string_vector(const std::string& Keyword) 
      const {
      return $self->value<std::vector<std::string> >(Keyword);
    }
    std::vector<double> value_double_vector(const std::string& Keyword) 
      const {
      return $self->value<std::vector<double> >(Keyword);
    }
  }
};
}

