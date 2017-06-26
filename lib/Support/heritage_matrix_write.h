#ifndef HERITAGE_MATRIX_WRITE_H
#define HERITAGE_MATRIX_WRITE_H
#include <string>
#include <blitz/array.h>
#include <map>

namespace FullPhysics {
  void heritage_matrix_write(const std::string& Fname, 
			     const std::string& File_id,
			     const blitz::Array<double, 2>& Arr,
			     const std::map<std::string, std::string>& 
			     Metadata = std::map<std::string, std::string>());
}
#endif
