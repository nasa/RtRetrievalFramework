#ifndef HERITAGE_FILE_H
#define HERITAGE_FILE_H
#include "printable.h"
#include "fp_exception.h"
#include "fp_time.h"
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <blitz/array.h>
#include <map>
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This class reads the heritage file formats. We read both the 
  configuration file and the matrix files (there are similar enough
  in format that it makes sense to combine these two).

  For the configuration files, we provide the values using a 
  "keyword path". By convention, this is just the keywords separated
  by a "/". For example, "CONTROL/input_file". Note that the keyword
  path is case insensitive, so "CONTROL/input_file" and 
  "control/INPUT_FILE" are the same.

  Some blocks are special in that there can be more than one of
  them. This includes "GAS", "INSTRUMENT", "AEROSOL" and "WINDOW". For
  blocks that we encounter more than one of, we access the data use a
  "index" number, which must be between 0 and the
  number_block(keyword). You can also view this data by adding the
  name to the keyword, so to find the moment file for the aerosol "IC"
  you can look at "PARAMETER_DEFINITION/AEROSOL/IC/moment_file". The
  two interfaces give the same data, use whichever is more
  convenient.

  For the various "HEADER" format (i.e., the matrix file), we just 
  use the keyword with no path, for example "Number_rows".

  We have an blitz Array "double". This is empty unless the file we
  read happens to contain matrix data. 

  We support conversion of the values read to a variety of
  formats. This is done by the "value" function. Currently supported
  conversions:

  - Any type that boost::lexical_cast<T> can convert a string to.
    In particular, double and int
  - A string. This includes stripping the quotes that may appear
    around the value
  - A boolean type. We translate the string "true" to true and "false"
    to false.
  - A std::vector<double>, which converts a list of doubles.
  - A std::vector<int>, which converts to a list of ints. This
    supports ranges in the entry, e.g. "1:4 7" which is returned as
    the list 1, 2, 3, 4, 7.
  - A std::vector<std::string>, which converts a set of strings where
    each is quoted.
  - A Time for a time stamp.

  Some of the values may be file names. The file and directory names
  are relative to the location of the heritage file. The routine
  file_value handles adding any necessary paths to the values found in
  the file.

  In addition, file_value will expand out environment variables like 
  "$(abscodir)".
*******************************************************************/

class HeritageFile : public Printable<HeritageFile> {
public:
  HeritageFile(const std::string& Fname);
  virtual ~HeritageFile() {}
  std::string file_value(const std::string& Keyword, int Block_index = 0) const;
  void parse_file(const std::string& Fname);
  void print(std::ostream& Os);

//-----------------------------------------------------------------------
/// Matrix data. This is an zero sized matrix if we aren't actually
/// reading a matrix file.
//-----------------------------------------------------------------------
  
  const blitz::Array<double, 2>& data() const {return data_;}

  int column_index(const std::string& Col_name) const;
  blitz::Array<double, 1> data(const std::string& Col_name) const;

//-----------------------------------------------------------------------
/// As an aid to error reporting, this adds the location that a
/// keyword was found to an Exception.
//-----------------------------------------------------------------------
  
  void add_location_to_exception(Exception& e,
				 const std::string& Keyword, 
				 int Block_index = 0) const
  {
    std::string k = Keyword;
    std::string v = get_value_or_exception(k, Block_index);
    e << "File: " << keyword_to_file.find(k)->second << "\n"
      << "Line: " << keyword_to_line.find(k)->second;
  }

//-----------------------------------------------------------------------
/// Return true if we found a value for the given Keyword path.
//-----------------------------------------------------------------------

  bool has_value(const std::string& Keyword) const
  {
    std::string k = Keyword;
    boost::to_lower(k);
    return (keyword_to_value.find(k) != keyword_to_value.end());
  }

//-----------------------------------------------------------------------
/// File name.
//-----------------------------------------------------------------------

  const std::string& file_name() const {return file_name_;}

//-----------------------------------------------------------------------
/// Directory base. This the directory that file_value is relative to.
//-----------------------------------------------------------------------

  const std::string& directory_base() const {return dirbase; }

//-----------------------------------------------------------------------
/// Return number of blocks for a particular keyword (e.g.,
/// "PARAMETER_DEFINITION/GAS". If you call this with a block 
/// not found in the file, this is not an error. We just return 0 in
/// that case.
//-----------------------------------------------------------------------
  int number_block(const std::string& Keyword) const {
    std::string k = Keyword;
    boost::to_lower(k);
    if(section_to_count.count(k))
      return section_to_count.find(k)->second;
    else
      return 0;
  }

  template<class T> T value(const std::string& Keyword, 
			    int block_index = 0) const;
private:
  std::string file_name_;
  // This maps keywords to a value. As a convention, nest keywords
  // with a "/", for example "output/log_file". We map everything to
  // lower case.
  std::map<std::string, std::string> keyword_to_value;

  // File and line we read value at. This can be used for error
  // reporting
  std::map<std::string, std::string> keyword_to_file;
  std::map<std::string, int> keyword_to_line;

  std::map<std::string, int> section_to_count;
  std::map<std::string, std::map<int, std::string> > section_and_index_to_name;

  // Data, may be empty.
  blitz::Array<double, 2> data_;
  void parse_file_block_count(const std::string& fname);
  std::string get_value_or_exception(std::string& k, 
				     int Block_index) const;

  // Directory files are relative to
  std::string dirbase;
};

//-----------------------------------------------------------------------
/// Return the value for the given Keyword path, cast to the given type T.
//-----------------------------------------------------------------------

template<class T> inline T HeritageFile::value(const 
		       std::string& Keyword, int Block_index) const
{
  std::string k = Keyword;
  std::string v = get_value_or_exception(k, Block_index);
  try {
    return boost::lexical_cast<T>(v);
  } catch(boost::bad_lexical_cast) {
    Exception e(
"Error reading file. The value for " + Keyword + " is not\n"
"the expected type\n");
    e << "File: " << keyword_to_file.find(k)->second << "\n"
      << "Line: " << keyword_to_line.find(k)->second;
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Return the value for the given Keyword path, cast to the given type T.
//-----------------------------------------------------------------------

template<> inline Time
HeritageFile::value<Time>(const 
		       std::string& Keyword, int Block_index) const
{
  std::string k = Keyword;
  std::string s = get_value_or_exception(k, Block_index);
  try {
    return Time::parse_time(s);
  } catch(boost::bad_lexical_cast) {
    Exception e(
"Error reading file. The value for " + Keyword + " is not\n"
"a datetime\n");
    e << "File: " << keyword_to_file.find(k)->second << "\n"
      << "Line: " << keyword_to_line.find(k)->second;
    throw e;
  }
}

//-----------------------------------------------------------------------
/// Return the value for the given Keyword path, cast to the given type T.
//-----------------------------------------------------------------------

template<> inline std::string HeritageFile::value<std::string>(const 
		       std::string& Keyword, int Block_index) const
{
  std::string k = Keyword;
  std::string v = get_value_or_exception(k, Block_index);
  boost::smatch m;
  if(boost::regex_match(v, m, boost::regex("\\s*\"?([^\"]*)\"?\\s*")))
    return m[1];
  else
    return "";
}

//-----------------------------------------------------------------------
/// Return the value for the given Keyword path, cast to the given type T.
//-----------------------------------------------------------------------

template<> inline 
std::vector<double> HeritageFile::value<std::vector<double> >
(const std::string& Keyword, int Block_index) const
{
  std::string k = Keyword;
  std::string v = get_value_or_exception(k, Block_index);
  std::vector<double> res;
  boost::sregex_token_iterator i(v.begin(), v.end(), boost::regex("\\s+"), -1);
  boost::sregex_token_iterator iend;
  if(i == iend || *i == "")		// Handle special case of an empty list
    return res;
  try {
    for(; i != iend; ++i)
      res.push_back(boost::lexical_cast<double>(*i));
  } catch(boost::bad_lexical_cast) {
    Exception e(
"Error reading file. The value for " + Keyword + " is not\n"
"the expected type\n");
    e << "File: " << keyword_to_file.find(k)->second << "\n"
      << "Line: " << keyword_to_line.find(k)->second;
    throw e;
  }
  return res;
}

//-----------------------------------------------------------------------
/// Return the value for the given Keyword path, cast to the given type T.
//-----------------------------------------------------------------------

template<> inline 
std::vector<int> HeritageFile::value<std::vector<int> >
(const std::string& Keyword, int Block_index) const
{
  std::string k = Keyword;
  std::string v = get_value_or_exception(k, Block_index);
  std::vector<int> res;
  boost::sregex_token_iterator i(v.begin(), v.end(), boost::regex("\\s+"), -1);
  boost::sregex_token_iterator iend;
  if(i == iend || *i == "")		// Handle special case of an empty list
    return res;
  try {
    boost::smatch m;
    for(; i != iend; ++i) {
      std::string ist = *i;
      if(boost::regex_match(ist,m, boost::regex("(\\d+):(\\d+)"))) {
	int rstart = boost::lexical_cast<int>(m[1]);
	int rend = boost::lexical_cast<int>(m[2]);
	for(int j = rstart; j <= rend; ++j)
	  res.push_back(j);
      } else if(boost::regex_match(ist,m, boost::regex("(\\d+)"))) {
	res.push_back(boost::lexical_cast<int>(m[1]));
      } else if (ist.length() > 0) {
	Exception e("Error reading file. The value for " + Keyword + " contains\n"
		    "an expected token: \"" + *i + "\"\n");
	e << "File: " << keyword_to_file.find(k)->second << "\n"
	  << "Line: " << keyword_to_line.find(k)->second;
	throw e;
      }
    }
  } catch(boost::bad_lexical_cast) {
    Exception e("Error reading file. The value for " + Keyword + " is not\n"
		"the expected type\n");
    e << "File: " << keyword_to_file.find(k)->second << "\n"
      << "Line: " << keyword_to_line.find(k)->second;
    throw e;
  }
  return res;
}

//-----------------------------------------------------------------------
/// Return the value for the given Keyword path, cast to the given type T.
//-----------------------------------------------------------------------

template<> inline 
std::vector<std::string> HeritageFile::value<std::vector<std::string> >
(const std::string& Keyword, int Block_index) const
{
  std::string k = Keyword;
  std::string v = get_value_or_exception(k, Block_index);
  boost::sregex_token_iterator icur;
  boost::sregex_token_iterator iend;

  // Try parsing strings in quotes seperated by spaces, ex:
  // "LABEL_1" "LABEL_2" "LABEL_3"
  icur = boost::sregex_token_iterator(v.begin(), v.end(), 
				      boost::regex("\\s*\"([^\"]*)\""), 
				      1);

  // If nothing was found from that parsing try just unquoted word blocks, ex:
  // LABEL_1 LABEL2 LABEL_3
  if (icur == iend)
    icur = boost::sregex_token_iterator(v.begin(), v.end(), 
					boost::regex("\\s*(\\w+)"), 
					1);
  if(icur == iend || *icur == "") // Handle special case of an empty list
    return std::vector<std::string>();
  return std::vector<std::string>(icur, iend);
}

//-----------------------------------------------------------------------
/// Return the value for the given Keyword path, cast to the given type T.
//-----------------------------------------------------------------------

template<> inline bool HeritageFile::value<bool>(const 
		 std::string& Keyword, int Block_index) const
{
  std::string k = Keyword;
  std::string t = get_value_or_exception(k, Block_index);
  boost::to_lower(t);
  if(t != "true" && t != "false") {
    Exception e(
"Error reading file. The value for " + Keyword + " is not\n"
"'true' or 'false'\n");
    e << "File: " << keyword_to_file.find(k)->second << "\n"
      << "Line: " << keyword_to_line.find(k)->second;
    throw e;
  }
  return(t == "true");
}
}
#endif
