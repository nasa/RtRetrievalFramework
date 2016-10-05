#ifndef HDF_FILE_GENERATING_H
#define HDF_FILE_GENERATING_H
#include "hdf_file.h"
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/****************************************************************//**
  To avoid creating files when an error occurs, we create the file
  with the name ".generating" appended. Only after the file has been
  fully rewritten to we rename it without the ".generating" extension
  added. 

  This class handles this bit of logic.
*******************************************************************/

class HdfFileGenerating : public Printable<HdfFileGenerating> {
public:
//-----------------------------------------------------------------------
/// Create a HDF file with the given name. This should *not* have the
/// ".generating" added, we add this in this class.
//-----------------------------------------------------------------------

  HdfFileGenerating(const std::string& Fname) : fname (Fname) {}
  virtual ~HdfFileGenerating() { close(); }
  void close();
  HdfFile& hdf_file();
  void abandon();
  void print(std::ostream& Os) const;
private:
  std::string fname;
  boost::shared_ptr<HdfFile> hdf_file_;
};
}
#endif
