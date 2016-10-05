#ifndef DIR_CHANGE_H
#define DIR_CHANGE_H
#include <string>
#include <unistd.h>
#include <fcntl.h>
#include "fp_exception.h"

namespace FullPhysics {

//-----------------------------------------------------------------------
/// Utility class. This changes to a new directory, and on destruction
/// changes back.
//-----------------------------------------------------------------------

class DirChange {
public:
  DirChange(const std::string& newdir)
  {
    dirhandle = open(".", O_RDONLY);
    int status = chdir(newdir.c_str());
    if(status != 0) {
      std::stringstream err_msg;
      err_msg << "Could not change to directory: " << newdir;
      throw Exception(err_msg.str());
    }
  }
  ~DirChange() 
  {
    int status = fchdir(dirhandle);
    close(dirhandle);
    if(status != 0)
      throw Exception("Call to fchdir failed");
  }
private:
  int dirhandle;
};
}
#endif
