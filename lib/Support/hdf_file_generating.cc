#include "hdf_file_generating.h"

using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Close a file, changing the name to the final name. This is
/// automatically called by the destructor, so you don't need to call
/// this directly.
//-----------------------------------------------------------------------

void HdfFileGenerating::close()
{
  if(hdf_file_) {
    hdf_file_->close();
    hdf_file_.reset();
    int status = rename((fname + ".generating").c_str(), fname.c_str());
    if(status != 0)
      throw Exception("Trouble renaming the file " + fname + 
		      ".generating to " + fname);
  }
}

HdfFile& HdfFileGenerating::hdf_file()
{
  if(!hdf_file_) {
  // Ok if this fails, we just want to remove any file that might be
  // in the way
    remove(fname.c_str());
    // Create with ".generating", which we remove when done creating file.
    hdf_file_.reset(new HdfFile(fname + ".generating", HdfFile::CREATE));
  }
  return *hdf_file_;
}

//-----------------------------------------------------------------------
/// Abandon a file and clean it up. This would normally be called if
/// an error condition occurs and we can't finish creating the file.
//-----------------------------------------------------------------------

void HdfFileGenerating::abandon()
{
  if(hdf_file_) {
    // Try to close first, but don't fail if closing fails.
    try {
      hdf_file_->close();
    } catch(...) {
    }
    hdf_file_.reset();
    // Don't bother checking return code, if removing fails we don't
    // have anything else to try.
    remove((fname + ".generating").c_str());
  }
}

void HdfFileGenerating::print(std::ostream& Os) const
{
  Os << "HdfFileGenerating: \n" 
     << "  File name: " << fname << "\n";
}
