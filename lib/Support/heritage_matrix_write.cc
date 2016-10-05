#include "heritage_matrix_write.h"
#include "fstream_compress.h"
#include <time.h>

//-----------------------------------------------------------------------
/// This writes a matrix out with the old Matrix IO format. This will
/// automatically compress if Fname ends with a ".gz". This can be
/// read by HeritageFile.
///
/// You can optionally pass in metadata, which is given as key/value
/// pairs. We just support values of strings, you need to convert
/// other values to a string representation before passing it. This is
/// a bit of a limited interface, but it is sufficient for what we
/// need to use this for.
//-----------------------------------------------------------------------

void FullPhysics::heritage_matrix_write(const std::string& Fname, 
			   const std::string& File_id,
			   const blitz::Array<double, 2>& Arr,
			   const std::map<std::string, std::string>& Metadata)
{
  using namespace blitz;
  time_t tm;
  time(&tm);			// Get time, to report in headers that
  OstreamCompress out(Fname);
  out.precision(9);
  out << std::scientific;
  out << "begin HEADER\n"
      << "  File_ID       = " << File_id << "\n"
      << "  File_Creation = " << ctime(&tm)
      << "  BuildId       = exported\n";
  if(Arr.rows() > 0)
    out << "  File_Type     = Matrix\n"
	<< "  Num_Rows      = " << Arr.rows() << "\n"
	<< "  Num_Columns   = " << Arr.cols() << "\n";
  for(std::map<std::string, std::string>::const_iterator i = Metadata.begin();
      i != Metadata.end(); ++i)
    out << "  " << i->first << " = " << i->second << "\n";
  out << "end HEADER\n";
  for(int i = 0; i < Arr.rows(); ++i) {
    for(int j = 0; j < Arr.cols(); ++j)
      out << Arr(i, j) << " ";
    out << "\n";
  }
}
