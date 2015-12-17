#include "heritage_file.h"
#include "fstream_compress.h"
#include "environment_substitute.h"
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <sys/stat.h>
#include "unistd.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
typedef const blitz::Array<double, 2>& (HeritageFile::*d1)(void) const;
typedef blitz::Array<double, 1> (HeritageFile::*d2)(const std::string&) const;

double value_double(const HeritageFile& f, const std::string keyword) {
  return f.value<double>(keyword);
}

int value_int(const HeritageFile& f, const std::string keyword) {
  return f.value<int>(keyword);
}

std::string value_string(const HeritageFile& f, const std::string keyword) {
  return f.value<std::string>(keyword);
}

#include "register_lua.h"
REGISTER_LUA_CLASS(HeritageFile)
.def(luabind::constructor<std::string>())
.def("data", ((d1) &HeritageFile::data))
.def("data", ((d2) &HeritageFile::data))
.def("column_index", &HeritageFile::column_index)
.def("has_value", &HeritageFile::has_value)
.def("value_double", &value_double)
.def("value_int", &value_int)
.def("value", &value_string)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Read the given heritage file.
//-----------------------------------------------------------------------

HeritageFile::HeritageFile(const std::string& Fname) 
:file_name_(Fname)
{
  // The various files in the configuration file are relative to 
  // the directory of the run file. So break into file name and
  // directory base.
  size_t t = Fname.find_last_of("/");
  if(t != std::string::npos)
    dirbase = Fname.substr(0, t);
  else
    dirbase = ".";
  dirbase += "/";
  parse_file(Fname);
}

//-----------------------------------------------------------------------
/// Make a pass through the data and just count the number of
/// blocks. We do this because we treat blocks with a count of 1
/// different than with a count > 1.
//-----------------------------------------------------------------------

void HeritageFile::parse_file_block_count(const std::string& fname)
{
  IstreamCompress ifn(fname);
  if(!ifn.good())
    throw Exception("Trouble reading file " + fname);
  ifn.exceptions(std::ifstream::badbit);
  std::string ln;
  // We handle "HEADER" blocks a little differently. We only have one
  // in a file, and after a HEADER block we may have matrix data. We
  // also leave HEADER outside of the keyword. So track when we enter
  // a header.
  bool have_header = false;
  bool done_with_header = false;
  std::vector<std::string> keyword;
  int line = 0;			// Keep track of the line we are
				// reading so we can report it in any
				// error message
  try {
    while(!ifn.eof() && ifn.good() && !done_with_header) {
      getline(ifn, ln);
      line++;
      if(ifn.good()) {

// -----------------------------------------------------------------------
// Strip off comments, and replace environment variables.
// -----------------------------------------------------------------------

	ln = ln.substr(0,ln.find("#"));
	ln = environment_substitute(ln);
	boost::smatch m;

// -----------------------------------------------------------------------
// Check if we are:
//   1) In a HEADER block.
// -----------------------------------------------------------------------

	if(boost::regex_match(ln, m, boost::regex("(?i) *begin +header *"))) {
	  have_header = true;
	  if(keyword.size() !=0)
	    throw Exception(
"Error reading file. Have a HEADER block nested inside of another block"
			    );

// -----------------------------------------------------------------------
//   2) Ending a HEADER block.
// -----------------------------------------------------------------------

	} else if(boost::regex_match(ln, m, boost::regex("(?i) *end +header *"))) {
	  if(!have_header)
	    throw Exception(
 "Error reading file. End of header seen before header started"
			    );
	  done_with_header = true;

// -----------------------------------------------------------------------
//   3) Beginning a section
// -----------------------------------------------------------------------

	} else if(boost::regex_match(ln, m, boost::regex("(?i) *begin +(\\w+) *"))) {
	  std::string t = m[1];
	  keyword.push_back(t);
	  std::string valkey = keyword[0];
	  for(int i = 1; i < (int) keyword.size(); ++i)
	    valkey = valkey + "/" + keyword[i];
	  boost::to_lower(valkey);
	  if(section_to_count.count(valkey) == 0)
	    section_to_count[valkey] = 1;
	  else 
	    section_to_count[valkey] += 1;

// -----------------------------------------------------------------------
//   4) Ending a section
// -----------------------------------------------------------------------

	} else if(boost::regex_match(ln, m, boost::regex("(?i) *end +(\\w+) *"))) {
	  std::string t = m[1];
	  if(keyword.size() < 1 ||
	     keyword.back() != t)
	    throw Exception(
"Error reading file. Expected to end section for\n"
"'" + keyword.back() + "' instead section end was for '" +  t + "'");
	  keyword.pop_back();

// -----------------------------------------------------------------------
//   5) Giving the name in a section
// -----------------------------------------------------------------------

	} else if(boost::regex_match(ln, m, boost::regex("(?i) *name *= *(.+?) *"))) {
	  if(!have_header && keyword.size() >= 1) {
	    std::string valkey = keyword[0];
	    for(int i = 1; i < (int) keyword.size(); ++i)
	      valkey = valkey + "/" + keyword[i];
	    boost::to_lower(valkey);
	    section_and_index_to_name[valkey][section_to_count[valkey] - 1] = 
	      m[1];
	  }

// -----------------------------------------------------------------------
//   6) Skip other data for now, we handle that in second pass.
// -----------------------------------------------------------------------

	} else {
	}
      }
    }
    if(keyword.size() != 0)
      throw Exception(
"Error reading file. Missing end of section for '" + 
keyword.back() + "'");
  } catch(Exception& exc) {
    exc << "\nFile: " << fname << "\n"
	<< "Line: " << line;
    throw;
  }
}

//-----------------------------------------------------------------------
/// Parse the given file, and add the keyword information from that
/// file to the keyword list.
//-----------------------------------------------------------------------

void HeritageFile::parse_file(const std::string& fname)
{
  parse_file_block_count(fname);
  IstreamCompress ifn(fname);
  if(!ifn.good())
    throw Exception("Trouble reading file " + fname);
  ifn.exceptions(std::ifstream::badbit);
  std::string ln;
  // We handle "HEADER" blocks a little differently. We only have one
  // in a file, and after a HEADER block we may have matrix data. We
  // also leave HEADER outside of the keyword. So track when we enter
  // a header.
  bool have_header = false;
  bool done_with_header = false;
  std::vector<std::string> keyword;
  std::map<std::string, int> cur_idx;
  int line = 0;			// Keep track of the line we are
				// reading so we can report it in any
				// error message
  try {
    while(!ifn.eof() && ifn.good() && !done_with_header) {
      getline(ifn, ln);
      line++;
      if(ifn.good()) {

// -----------------------------------------------------------------------
// Strip off comments, and replace environment variables.
// -----------------------------------------------------------------------

	ln = ln.substr(0,ln.find("#"));
	ln = environment_substitute(ln);
	boost::smatch m;

// -----------------------------------------------------------------------
// Check if we are:
//   1) In a HEADER block.
// -----------------------------------------------------------------------

	if(boost::regex_match(ln, m, boost::regex("(?i) *begin +header *"))) {
	  have_header = true;
	  if(keyword.size() !=0)
	    throw Exception(
"Error reading file. Have a HEADER block nested inside of another block"
			    );

// -----------------------------------------------------------------------
//   2) Ending a HEADER block.
// -----------------------------------------------------------------------

	} else if(boost::regex_match(ln, m, boost::regex("(?i) *end +header *"))) {
	  if(!have_header)
	    throw Exception(
 "Error reading file. End of header seen before header started"
			    );
	  done_with_header = true;

// -----------------------------------------------------------------------
//   3) Beginning a section
// -----------------------------------------------------------------------

	} else if(boost::regex_match(ln, m, boost::regex("(?i) *begin +(\\w+) *"))) {
	  std::string t = m[1];
	  keyword.push_back(t);
	  std::string valkey = keyword[0];
	  for(int i = 1; i < (int) keyword.size(); ++i)
	    valkey = valkey + "/" + keyword[i];
	  boost::to_lower(valkey);
	  if(cur_idx.count(valkey) == 0)
	    cur_idx[valkey] = 0;
	  else 
	    cur_idx[valkey] += 1;
	  

// -----------------------------------------------------------------------
//   4) Ending a section
// -----------------------------------------------------------------------

	} else if(boost::regex_match(ln, m, boost::regex("(?i) *end +(\\w+) *"))) {
	  std::string t = m[1];
	  if(keyword.size() < 1 ||
	     keyword.back() != t)
	    throw Exception(
"Error reading file. Expected to end section for\n"
"'" + keyword.back() + "' instead section end was for '" +  t + "'");
	  keyword.pop_back();

// -----------------------------------------------------------------------
//   5) Setting a value. We allow this to be empty, for keywords set
//   to null.
// -----------------------------------------------------------------------

	} else if(boost::regex_match(ln, m, boost::regex(" *(.+?) *= *(.*?) *"))) {
	  if(have_header) {
	    std::string valkey = m[1];
	    boost::to_lower(valkey);
	    keyword_to_value[valkey] = m[2];
	    keyword_to_file[valkey] = fname;
	    keyword_to_line[valkey] = line;
	  } else {
	    if(keyword.size() < 1)
	      throw Exception(
"Error reading file. Value assigned outside of a\n"
"keyword section");

	    // Set up the values using the keyword path + index number
	    // (if we have more than one block)
	    std::string valkey = keyword[0];
	    for(int i = 1; i < (int) keyword.size(); ++i)
	      valkey = valkey + "/" + keyword[i];
	    boost::to_lower(valkey);
	    if(section_to_count.count(valkey) > 0 &&
	       section_to_count[valkey] > 1)
	      valkey += "/" + 
		boost::lexical_cast<std::string>(cur_idx[valkey]);
	    valkey = valkey + "/" + m[1];
	    boost::to_lower(valkey);
	    keyword_to_value[valkey] = m[2];
	    keyword_to_file[valkey] = fname;
	    keyword_to_line[valkey] = line;

	    // Record exactly the same information, but have keyword
	    // path + "name" of the block. This may be a more
	    // convenient way for the customer of this class to access
	    // the data.

	    valkey = keyword[0];
	    for(int i = 1; i < (int) keyword.size(); ++i)
	      valkey = valkey + "/" + keyword[i];
	    boost::to_lower(valkey);
	    if(section_and_index_to_name.count(valkey) > 0 &&
	       section_and_index_to_name[valkey].count(cur_idx[valkey]) > 0) {
	      valkey = valkey + "/" + 
		section_and_index_to_name[valkey][cur_idx[valkey]];
	      valkey = valkey + "/" + m[1];
	      boost::to_lower(valkey);
	      keyword_to_value[valkey] = m[2];
	      keyword_to_file[valkey] = fname;
	      keyword_to_line[valkey] = line;
	    }
	  }

// -----------------------------------------------------------------------
//   6) Have a stray line.
// -----------------------------------------------------------------------

        } else if(boost::regex_search(ln, m, boost::regex("\\w"))) {
	  throw Exception("Error reading file. Unrecognized line");
	}
      }
    }
    if(keyword.size() != 0)
      throw Exception("Error reading file. Missing end of section for '" + 
		      keyword.back() + "'");

// -----------------------------------------------------------------------
/// Read matrix data, if we have any.
// -----------------------------------------------------------------------

    if(have_header && has_value("num_rows")) {
      int nr = value<int>("num_rows");
      int nc = value<int>("num_columns");
      data_.resize(nr, nc);
      int i = 0;
      while(!ifn.eof() && ifn.good() && i < nr) {
	getline(ifn, ln);
	line++;
	if(ifn.good()) {

// -----------------------------------------------------------------------
// Strip off comments, and replace environment variables.
// -----------------------------------------------------------------------

	  ln = ln.substr(0,ln.find("#"));
	  ln = environment_substitute(ln);

// -----------------------------------------------------------------------
// Replace any "d" with "e", this is the difference in the way values
// are expressed in fortran vs. C.
// -----------------------------------------------------------------------
	  
	  boost::replace_all(ln, std::string("d"), std::string("e"));
	  boost::smatch m;
	  if(boost::regex_search(ln, m, boost::regex("\\w"))) {
	    std::istringstream is(ln);
	    is.exceptions(std::ifstream::badbit);
	    for(int j = 0; j < nc; ++j)
	      is >> data_(i, j);
	    if(is.fail())
	      throw Exception(
"Error reading file. Trouble processing a row of data");
	    ++i;
	  }
	}
      }
      if(i < nr)
	throw Exception("Error reading file. Not enough lines of data");
    }
  } catch(Exception& exc) {
    exc << "\nFile: " << fname << "\n"
	<< "Line: " << line;
    throw;
  }
}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

void HeritageFile::print(std::ostream& Os)
{
  Os << "Heritage file:\n"
     << "  File name: " << file_name_;
}


//-----------------------------------------------------------------------
/// Get a value, or throw an exception if we don't find a value for 
/// the keyword. We may update the keyword in place to include the
/// Block_index. 
//-----------------------------------------------------------------------

std::string HeritageFile::get_value_or_exception(std::string& k, 
						 int Block_index) const
{
  std::string korig = k;	// Use for error reporting.
  boost::to_lower(k);
  // Find enclosing block, if any. If we have an enclosing block and
  // more than one of them, then add the block index to the keyword
  // string. 
  size_t i = k.find_last_of('/');
  if(i != std::string::npos) {
    std::string s1 = k.substr(0, i);
    std::string s2 = k.substr(i + 1);
    if(number_block(s1) > 0){	// If 0, drop through. We catch the
				// error after this block.
      range_check(Block_index, 0, number_block(s1));
      if(number_block(s1) > 1)
	k = s1 + "/" + boost::lexical_cast<std::string>(Block_index) + "/" +
	  s2;
    }
  }
  if(!has_value(k))
    throw Exception(
"Error reading file. Expected to find keyword " + korig + "\n"
"but didn't find in file.\n"
"File:   " + file_name());
  return keyword_to_value.find(k)->second;
}

//-----------------------------------------------------------------------
/// Return a file or directory for the given keyword path. This
/// handles converting the file or directory relative to the directory
/// that the run file is in.
///
/// We also expand environment variables given like "$(abscodir)"
//-----------------------------------------------------------------------

std::string HeritageFile::file_value(const std::string& Keyword, 
				     int Block_index) 
  const 
{ 
  std::string f = value<std::string>(Keyword, Block_index);
  
  std::string newfn;
  if(f[0] == '/')	   // Just return f if it is an absolute path.
    newfn = f;
  else			   // Otherwise, add base path.
    newfn = dirbase + f;
  
  struct stat buf;
  if (stat(newfn.c_str(), &buf) != 0) {
    Exception errmsg;
    errmsg << "File referenced by keyword: " << Keyword << ", "
	   << "for block index: " << Block_index << " does not exist as: "
	   << newfn;
    throw errmsg;
  }
  
  return newfn;
}

//-----------------------------------------------------------------------
/// For files with a Label attribute, this is the index of the named 
/// column. This relies upon having a labels keyword in the file.
//-----------------------------------------------------------------------

int HeritageFile::column_index(const std::string& Col_name) const {
  std::string colnamelc = Col_name;
  boost::to_lower(colnamelc);
  std::vector<std::string> lb = value<std::vector<std::string> >("labels");
  BOOST_FOREACH(std::string& t, lb)
    boost::to_lower(t);
  std::vector<std::string>::const_iterator i = 
    std::find(lb.begin(), lb.end(), colnamelc);
  if(i ==lb.end())
    return -1;
  else
    return i - lb.begin();
}

//-----------------------------------------------------------------------
/// For files with a Label attribute, this is the column that has the
/// given label.
//-----------------------------------------------------------------------

blitz::Array<double, 1> HeritageFile::data(const std::string& Col_name) 
  const
{
  int ind = column_index(Col_name);
  if(ind < 0)
    throw Exception("Did not find column '" + Col_name + "' in the file " + 
		    file_name());
  return data()(blitz::Range::all(), ind);
}
