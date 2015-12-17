#ifndef HDF_FILE_H
#define HDF_FILE_H
#include "printable.h"
#include "fp_exception.h"
#include "linear_algebra.h"
#include "array_with_unit.h"
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>
#include <H5Cpp.h>
#include <vector>
#include <stdint.h>
#include <string.h>

namespace FullPhysics {
/****************************************************************//**
  This class reads and writes a HDF5 file. Note that this is just a
  thin layer on top of the HDF 5 libraries to make the file operations
  we need to do in Level 2 Full Physics easier. There are many other
  things that can be done with a HDF 5 than what this class exposes.

  Note that because it is what is used by Level 2 product, we produce
  data in 32 bit (either 32 bit integer or 32 bit floating point). On
  a 64 bit system, the underlying double and int are larger. We map
  between these types as needed transparently.

  HDF supports both fixed length strings and variable length
  strings. We have the need to write both variable length strings and
  fixed length strings. As a convention, if std::string are passed we
  write as variable length. If const char* is passed, we determine the
  fixed length needed to accommodate the largest string passed and
  write the data as fixed strings of that length (padding shorter
  strings with spaces). This is an arbitrary decision, but it allows
  us to write both types. 

  Note that in what is a fairly odd convention, we add a trailing '\0'
  in our fixed length string, so there is one extra character.
*******************************************************************/

class HdfFile : public Printable<HdfFile> {
public:
  enum Mode {READ, CREATE, READ_WRITE};
  HdfFile(const std::string& Fname, Mode M = READ);
  virtual ~HdfFile() {}

//-----------------------------------------------------------------------
/// Check to see if an object (such as a Dataset) is in the file.
//-----------------------------------------------------------------------
  bool has_object(const std::string& Objname) const
  { return is_present(Objname, *h); }

//-----------------------------------------------------------------------
/// Check to see if a attribute is in the file.
//-----------------------------------------------------------------------
  bool has_attribute(const std::string& Aname) const;

  Unit read_units(const std::string& Dataname) const;
  template<int D> blitz::TinyVector<int,D> read_shape(
         const std::string& Dataname) const;
  template<class T, int D> blitz::Array<T, D> read_field(
	 const std::string& Dataname) const;
  template<class T, int D> ArrayWithUnit<T, D> read_field_with_unit(
	 const std::string& Dataname) const;
  template<class T, int D> ArrayWithUnit<T, D> read_field_with_unit(
        const std::string& Dataname, const Unit& Default_unit) const;
  template<class T, int D> ArrayWithUnit<T, D> read_field_with_unit(
	const std::string& Dataname,
	const blitz::TinyVector<int,D>& Start, 
	const blitz::TinyVector<int,D>& Size) const;
  template<class T, int D> ArrayWithUnit<T, D> read_field_with_unit(
        const std::string& Dataname, const Unit& Default_unit,
	const blitz::TinyVector<int,D>& Start, 
	const blitz::TinyVector<int,D>& Size) const;
  template<class T, int D> blitz::Array<T, D> read_field(
	 const std::string& Dataname, 
	 const blitz::TinyVector<int,D>& Start, 
	 const blitz::TinyVector<int,D>& Size) const;
  template<class T> T read_field(
	 const std::string& Dataname) const;
  template<class T, int D> blitz::Array<T, D> 
  read_attribute(const std::string& Aname) const;
  template<class T> T read_attribute(const std::string& Aname) const;
  template<class T, int D> void write_attribute(const std::string& Aname,
				const blitz::Array<T, D>& Data);
  template<class T> void write_attribute(const std::string& Aname,
				const T& Data);
  void write_attribute(const std::string& Aname, const char* Data);
  template<class T, int D> void write_field(const std::string& Dataname,
					    const blitz::Array<T, D>& Data);
  template<int D> inline void write_field(const std::string& Dataname, 
			  const blitz::Array<std::string, D>& Data);
  template<int D> inline void write_field(const std::string& Dataname, 
			  const blitz::Array<const char *, D>& Data);
  template<class T> void write_field(const std::string& Dataname, 
				     const T& Data);
  void write_field(const std::string& Dataname, 
		   const std::string& Data);
  void write_field(const std::string& Dataname, 
		   const char* Data);
  template<class T> void write_field(const std::string& Dataname, 
				     const T& Data,
				     H5::DataType P);
  template<class T, int D> void write_field(const std::string& Dataname,
					    const blitz::Array<T, D>& Data,
					    H5::DataType P);
  void dimension_metadata(const std::string& Name, 
			  const std::string& Description,
			  int Size);
  void shape_metadata(const std::string& Name, const std::string& Dim1);
  void shape_metadata(const std::string& Name, const std::string& Dim1,
		      const std::string& Dim2);
  void shape_metadata(const std::string& Name, const std::string& Dim1,
		      const std::string& Dim2, const std::string& Dim3);

//-----------------------------------------------------------------------
/// Mode file was opened with.
//-----------------------------------------------------------------------

  Mode mode() const {return mode_;}

//-----------------------------------------------------------------------
/// Close the underlying file. This is automatically done by the
/// destructor, so you only need to call this if you want to force a
/// close (e.g., for a unit test)
//-----------------------------------------------------------------------

  void close() { h->close(); }

//-----------------------------------------------------------------------
/// File name
//-----------------------------------------------------------------------
  const std::string& file_name() const {return fname;}

//-----------------------------------------------------------------------
/// Return true if the given file is an HDF file.
//-----------------------------------------------------------------------

  static bool is_hdf(const std::string& Fname) 
  { 
    try {
      return H5::H5File::isHdf5(Fname); 
    } catch(const H5::FileIException& e) {
      FullPhysics::Exception err;
      err << "Error determining if file: "
	  << Fname << " is an HDF file: "
	  << e.getDetailMsg();
      throw err;
    }
  }

  // Map native type to type of blitz array
  template<class T> H5::PredType pred_arr() const;
  // Map type of blitz array to type of data in file.
  template<class T> H5::PredType pred_data() const;
  void print(std::ostream& Os) const;
  H5::H5File& h5_file() { return *h; };
  const H5::H5File& h5_file() const { return *h; };
private:
  boost::shared_ptr<H5::H5File> h;
  std::string fname;
  Mode mode_;
  void create_group_if_needed(const std::string& Dataname, 
			      H5::CommonFG& Parent);
  H5::Attribute open_attribute(const std::string& Aname) const;
  H5::Attribute create_attribute(const std::string& Aname, 
	 const H5::DataSpace& Ds, const H5::DataType& P);
  bool is_group(const std::string& Objname) const;
  bool is_present(const std::string& Objname, 
		  const H5::CommonFG& Parent) const;
  void write_type(const std::string& Dataname, const H5::DataType& P);
};

template<> inline H5::PredType HdfFile::pred_arr<int>() const 
{return H5::PredType::NATIVE_INT;}
template<> inline H5::PredType HdfFile::pred_arr<int64_t>() const 
{return H5::PredType::NATIVE_INT64;}
template<> inline H5::PredType HdfFile::pred_arr<double>() const 
{return H5::PredType::NATIVE_DOUBLE;}
template<> inline H5::PredType HdfFile::pred_arr<float>() const 
{return H5::PredType::NATIVE_FLOAT;}
template<> inline H5::PredType HdfFile::pred_data<int>() const 
{return H5::PredType::NATIVE_INT32;}
template<> inline H5::PredType HdfFile::pred_data<int64_t>() const 
{return H5::PredType::NATIVE_INT64;}
template<> inline H5::PredType HdfFile::pred_data<double>() const
{return H5::PredType::NATIVE_DOUBLE;}
template<> inline H5::PredType HdfFile::pred_data<float>() const
{return H5::PredType::NATIVE_FLOAT;}

//-----------------------------------------------------------------------
/// Read the given attribute attached to a group or dataset.
//-----------------------------------------------------------------------

template<class T, int D> inline blitz::Array<T, D> 
HdfFile::read_attribute(const std::string& Aname) const
{
  try {
  using namespace H5;
  Attribute a = open_attribute(Aname);
  DataSpace ds = a.getSpace();
  if(ds.getSimpleExtentNdims() != D) {
    Exception e;
    e << "Attribute " << Aname << " does not have the expected rank of " 
      << D;
    throw e;
  }
  hsize_t dims[D];
  ds.getSimpleExtentDims(dims, NULL);
  blitz::TinyVector<int,D> dims2;
  for(int i = 0; i < D; ++i)
    dims2(i) = (int) dims[i];
  blitz::Array<T, D> res(dims2);
  a.read(pred_arr<T>(), res.dataFirst());
  return res;
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While reading attribute " << Aname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

//-----------------------------------------------------------------------
/// Read the given attribute attached to a group or dataset.
//-----------------------------------------------------------------------

template<class T> inline T 
HdfFile::read_attribute(const std::string& Aname) const
{
  return read_attribute<T, 1>(Aname)(0);
}

//-----------------------------------------------------------------------
/// Read the given attribute attached to a group or dataset.
//-----------------------------------------------------------------------

template<> inline std::string HdfFile::read_attribute<std::string>(const 
	    std::string& Aname) const
{
  try {
  using namespace H5;
  Attribute a = open_attribute(Aname);
  std::string res;
  a.read(a.getStrType(),res);
  // Fixed length strings may have trailing spaces. Truncate them.
  size_t t = res.find_last_not_of(" ");
  if(t != std::string::npos)
    res.erase(t + 1);
  return res;
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While reading attribute " << Aname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

//-----------------------------------------------------------------------
/// Read the given attribute attached to a group or dataset.
//-----------------------------------------------------------------------

template<> inline std::vector<std::string> 
HdfFile::read_attribute<std::vector<std::string> >(const 
						   std::string& Aname) const
{
  try {
  using namespace H5;
  Attribute a = open_attribute(Aname);
  DataSpace ds = a.getSpace();
  if(ds.getSimpleExtentNdims() != 1) {
    Exception e;
    e << "Attribute " << Aname << " does not have the expected rank of 1";
    throw e;
  }
  hsize_t dims[1];
  ds.getSimpleExtentDims(dims, NULL);
  int nelem = (int) dims[0];
  std::vector<char*> ptr(nelem);
  a.read(a.getStrType(),&(*ptr.begin()));
  std::vector<std::string> res(nelem);
  for(int i = 0; i < nelem; ++i) {
    res[i] = std::string(ptr[i]);
    // Fixed length strings may have trailing spaces. Truncate them.
    size_t t = res[i].find_last_not_of(" ");
    if(t != std::string::npos)
      res[i].erase(t + 1);
  }
  return res;
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While reading attribute " << Aname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

//-----------------------------------------------------------------------
// Somewhat annoyingly, C++ doesn't allow function partial
// specialization. Introduce a class as a roundabout way to do this,
// since we can do partial template specialization of a class.
//-----------------------------------------------------------------------
// Don't have Doxygen document this class.
/// @cond
template<int D> class HdfFilePartialDimHelper {
public:

blitz::TinyVector<int,D> read_shape(const HdfFile& hf,
				    const H5::H5File& h, 
				    const std::string& Dataname) const
{
  try {
    using namespace H5;
    DataSet d = h.openDataSet(Dataname);
    DataSpace ds = d.getSpace();
    if(ds.getSimpleExtentNdims() != D) {
      Exception e;
      e << "Dataset " << Dataname << " does not have the expected rank of " 
	<< D;
      throw e;
    }
    hsize_t dims[D];
    ds.getSimpleExtentDims(dims, NULL);
    blitz::TinyVector<int,D> dims2;
    for(int i = 0; i < D; ++i)
      dims2(i) = (int) dims[i];
    return dims2;
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While reading shape " << Dataname
       << " for the file '" << hf.file_name()
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

};

template<class T, int D> class HdfFilePartialSpecialHelper {
public:

blitz::Array<T, D> read_field(const HdfFile& hf,
			      const H5::H5File& h, 
			      const std::string& Dataname) const
{
  try {
  using namespace H5;
  DataSet d = h.openDataSet(Dataname);
  DataSpace ds = d.getSpace();
  if(ds.getSimpleExtentNdims() != D) {
    Exception e;
    e << "Dataset " << Dataname << " does not have the expected rank of " 
      << D;
    throw e;
  }
  hsize_t dims[D];
  ds.getSimpleExtentDims(dims, NULL);
  blitz::TinyVector<int,D> dims2;
  for(int i = 0; i < D; ++i)
    dims2(i) = (int) dims[i];
  blitz::Array<T, D> res(dims2);
  d.read(res.dataFirst(), hf.pred_arr<T>());
  return res;
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While reading field " << Dataname
       << " for the file '" << hf.file_name()
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

blitz::Array<T, D> read_field(const HdfFile& hf,
			      const H5::H5File& h, 
			      const std::string& Dataname,
			      const blitz::TinyVector<hsize_t,D>& Start, 
			      const blitz::TinyVector<hsize_t,D>& Size) const
{
  try {
  using namespace H5;
  DataSet d = h.openDataSet(Dataname);
  DataSpace ds = d.getSpace();
  if(ds.getSimpleExtentNdims() != D) {
    Exception e;
    e << "Dataset " << Dataname << " does not have the expected rank of " 
      << D;
    throw e;
  }
  DataSpace ms(D, &Size[0]);
  ds.selectHyperslab(H5S_SELECT_SET, &Size[0], &Start[0]);
  blitz::Array<T, D> res(Size);
  d.read(res.dataFirst(), hf.pred_arr<T>(), ms, ds);
  return res;
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While reading field " << Dataname
       << " for the file '" << hf.file_name()
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

};

template<int D> class HdfFilePartialSpecialHelper<std::string, D> {
public:
blitz::Array<std::string, D> read_field(const HdfFile& hf,
					const H5::H5File& h, 
					const std::string& Dataname) const
{
  try {
  using namespace H5;
  DataSet d = h.openDataSet(Dataname);
  DataSpace ds = d.getSpace();
  if(ds.getSimpleExtentNdims() != D) {
    Exception e;
    e << "Dataset " << Dataname << " does not have the expected rank of " 
      << D;
    throw e;
  }
  hsize_t dims[D];
  ds.getSimpleExtentDims(dims, NULL);
  blitz::TinyVector<int,D> dims2;
  for(int i = 0; i < D; ++i)
    dims2(i) = (int) dims[i];

  DataType dt = d.getDataType();
  blitz::Array<std::string, D> result_data(dims2);
  if(dt.isVariableStr()) {
    blitz::Array<const char*, D> read_data(result_data.shape());
    StrType st(PredType::C_S1, H5T_VARIABLE);
    d.read(read_data.dataFirst(), st);
    
    typename blitz::Array<const char*, D>::const_iterator i1;
    typename blitz::Array<std::string, D>::iterator i2 = result_data.begin();
    for(i1 = read_data.begin(); i1 != read_data.end(); ++i1, ++i2)
      *i2 = std::string(*i1);
  } else {
    int flat_size = 1;
    for(int i = 0; i < D; i++)
      flat_size *= dims2[i];
    flat_size *= dt.getSize();
    blitz::Array<char, 1> read_data(flat_size);
    d.read(read_data.dataFirst(), d.getStrType());

    char* i1 = read_data.dataFirst();
    typename blitz::Array<std::string, D>::iterator i2;
    for(i2 = result_data.begin(); i2 != result_data.end(); 
	++i2, i1 += dt.getSize()) {
      // By our fairly odd convention, this string may contain a
      // trailing '\0'. So we strip that off if we find it.
      std::string s(i1, i1 + dt.getSize());
      *i2 = std::string(s.c_str());
    }
  }
  return result_data;
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While reading field " << Dataname
       << " for the file '" << hf.file_name()
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

blitz::Array<std::string, D> read_field(const HdfFile& hf,
			const H5::H5File& h, 
			const std::string& Dataname,
		        const blitz::TinyVector<hsize_t,D>& Start, 
 		        const blitz::TinyVector<hsize_t,D>& Size) const
{
  try {
  using namespace H5;
  DataSet d = h.openDataSet(Dataname);
  DataSpace ds = d.getSpace();
  if(ds.getSimpleExtentNdims() != D) {
    Exception e;
    e << "Dataset " << Dataname << " does not have the expected rank of " 
      << D;
    throw e;
  }
  DataSpace ms(D, &Size[0]);
  ds.selectHyperslab(H5S_SELECT_SET, &Size[0], &Start[0]);
  DataType dt = d.getDataType();
  blitz::Array<std::string, D> result_data(Size);
  if(dt.isVariableStr()) {
    blitz::Array<const char*, D> read_data(result_data.shape());
    StrType st(PredType::C_S1, H5T_VARIABLE);
    d.read(read_data.dataFirst(), st, ms, ds);
    
    typename blitz::Array<const char*, D>::const_iterator i1;
    typename blitz::Array<std::string, D>::iterator i2 = result_data.begin();
    for(i1 = read_data.begin(); i1 != read_data.end(); ++i1, ++i2)
      *i2 = std::string(*i1);
  } else {
    int flat_size = 1;
    for(int i = 0; i < D; i++)
      flat_size *= Size(i);
    flat_size *= dt.getSize();
    blitz::Array<char, 1> read_data(flat_size);
    d.read(read_data.dataFirst(), d.getStrType(), ms, ds);

    char* i1 = read_data.dataFirst();
    typename blitz::Array<std::string, D>::iterator i2;
    for(i2 = result_data.begin(); i2 != result_data.end(); 
	++i2, i1 += dt.getSize()) {
      // By our fairly odd convention, this string may contain a
      // trailing '\0'. So we strip that off if we find it.
      std::string s(i1, i1 + dt.getSize());
      *i2 = std::string(s.c_str());
    }
  }
  return result_data;
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While reading field " << Dataname
       << " for the file '" << hf.file_name()
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

};

/// @endcond
  

//-----------------------------------------------------------------------
/// Reads the shape of a dataset
//-----------------------------------------------------------------------
template<int D> inline blitz::TinyVector<int,D>
HdfFile::read_shape(const std::string& Dataname) const
{
  HdfFilePartialDimHelper<D> helper;
  return helper.read_shape(*this, *h, Dataname);
}

//-----------------------------------------------------------------------
/// Read a given field.
//-----------------------------------------------------------------------

template<class T, int D> inline blitz::Array<T, D> 
HdfFile::read_field(const std::string& Dataname) const
{
  HdfFilePartialSpecialHelper<T,D> helper;
  return helper.read_field(*this, *h, Dataname);
}

//-----------------------------------------------------------------------
/// Read a given field, along with metadata describing the units.
//-----------------------------------------------------------------------

template<class T, int D> inline ArrayWithUnit<T, D> 
HdfFile::read_field_with_unit(const std::string& Dataname) const
{
  ArrayWithUnit<T, D> res;
  res.value.reference(read_field<T, D>(Dataname));
  res.units = read_units(Dataname);
  return res;
}

//-----------------------------------------------------------------------
/// Read a given field, along with metadata describing the units.
//  With a start and size given.
//-----------------------------------------------------------------------

template<class T, int D> inline ArrayWithUnit<T, D> 
HdfFile::read_field_with_unit(const std::string& Dataname, 
			      const blitz::TinyVector<int,D>& Start, 
			      const blitz::TinyVector<int,D>& Size) const
{
  ArrayWithUnit<T, D> res;
  res.data.reference(read_field<T, D>(Dataname, Start, Size));
  res.units = read_units(Dataname);
  return res;
}

//-----------------------------------------------------------------------
/// Read a given field, along with metadata describing the units.
/// Some data might be missing the units metadata field. For this
/// variant of read_field_with_unit, this won't be treated as an
/// error, but rather we will use the given default units instead.
//-----------------------------------------------------------------------

template<class T, int D> inline ArrayWithUnit<T, D> 
HdfFile::read_field_with_unit
(const std::string& Dataname, const Unit& Default_unit) const
{
  ArrayWithUnit<T, D> res;
  res.value.reference(read_field<T, D>(Dataname));
  if(!has_attribute(Dataname + "/Unit") &&
     !has_attribute(Dataname + "/Units")) 
    res.units = Default_unit;
  else
    res.units = read_units(Dataname);
  return res;
}

//-----------------------------------------------------------------------
/// Adds ability to specify Start and Size along with default units
//-----------------------------------------------------------------------

template<class T, int D> inline ArrayWithUnit<T, D> 
HdfFile::read_field_with_unit
(const std::string& Dataname, const Unit& Default_unit,
 const blitz::TinyVector<int,D>& Start, 
 const blitz::TinyVector<int,D>& Size) const
{
  ArrayWithUnit<T, D> res;
  res.value.reference(read_field<T, D>(Dataname, Start, Size));
  if(!has_attribute(Dataname + "/Unit") &&
     !has_attribute(Dataname + "/Units")) 
    res.units = Default_unit;
  else
    res.units = read_units(Dataname);
  return res;
}

//-----------------------------------------------------------------------
/// Read a given field. This reads a subset of the data, given by
/// Start and Size.
//-----------------------------------------------------------------------

template<class T, int D> blitz::Array<T, D> 
HdfFile::read_field(
	 const std::string& Dataname, 
	 const blitz::TinyVector<int,D>& Start, 
	 const blitz::TinyVector<int,D>& Size) const
{
  HdfFilePartialSpecialHelper<T,D> helper;
  blitz::TinyVector<hsize_t, D> start2, size2;
  start2 = Start;
  size2 = Size;
  return helper.read_field(*this, *h, Dataname, start2, size2);
}

//-----------------------------------------------------------------------
/// Read a given scalar field. This is a array of size 1, which we
/// return just the value in.
//-----------------------------------------------------------------------

template<class T> inline 
T HdfFile::read_field(const std::string& Dataname) const
{
  try {
    using namespace H5;
    DataSet d = h->openDataSet(Dataname);
    DataSpace ds = d.getSpace();
    T res;
    d.read(&res, pred_arr<T>());
    return res;
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While reading field " << Dataname
       << " for the file '" << file_name()
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

template<> inline 
std::string HdfFile::read_field(const std::string& Dataname) const
{
  return read_field<std::string, 1>(Dataname)(0);
}

//-----------------------------------------------------------------------
/// Write attribute to file
//-----------------------------------------------------------------------

template<class T, int D> inline
void HdfFile::write_attribute(const std::string& Aname,
			      const blitz::Array<T, D>& Data)
{
  try {
  using namespace H5;
  blitz::Array<T, D> data2 = to_c_order_const(Data);
  hsize_t dim[D];
  for(int i = 0; i < D; ++i)
    dim[i] = Data.extent(i);
  DataSpace ds(D, dim);
  Attribute a = create_attribute(Aname, ds, pred_data<T>());
  a.write(pred_arr<T>(), data2.dataFirst());
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While writing attribute " << Aname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

//-----------------------------------------------------------------------
/// Write attribute to file
//-----------------------------------------------------------------------

template<class T> inline 
void HdfFile::write_attribute(const std::string& Aname,
			      const T& Data)
{
  try {
  using namespace H5;
  hsize_t dim[1] = {1};
  DataSpace ds(1, dim);
  Attribute a = create_attribute(Aname, ds, pred_data<T>());
  a.write(pred_arr<T>(), &Data);
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While writing attribute " << Aname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

//-----------------------------------------------------------------------
/// Write attribute to file
//-----------------------------------------------------------------------

template<> inline 
void HdfFile::write_attribute<std::string>(const std::string& Aname,
			      const std::string& Data)
{
  try {
  using namespace H5;
  hsize_t dim[1] = {1};
  DataSpace ds(1, dim);
  StrType st(PredType::C_S1, H5T_VARIABLE);
  Attribute a = create_attribute(Aname, ds, st);
  a.write(st, Data);
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While writing attribute " << Aname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

//-----------------------------------------------------------------------
/// Write attribute to file
//-----------------------------------------------------------------------

template<> inline 
void HdfFile::write_attribute<std::vector<std::string> >(const 
    std::string& Aname, const std::vector<std::string>& Data)
{
  try {
  using namespace H5;
  int nelem = (int) Data.size();
  hsize_t dim[1] = {Data.size()};
  DataSpace ds(1, dim);
  StrType st(PredType::C_S1, H5T_VARIABLE);
  Attribute a = create_attribute(Aname, ds, st);
  std::vector<const char*> ptr(nelem);
  for(int i = 0; i < nelem; ++i)
    ptr[i] = Data[i].c_str();
  a.write(st, &(*ptr.begin()));
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While writing attribute " << Aname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

//-----------------------------------------------------------------------
/// Write a given field. If you give a group (e.g., "Group1/Field1"),
/// then we automatically create the group if it doesn't already exist.
//-----------------------------------------------------------------------

template<class T, int D> inline void
HdfFile::write_field(const std::string& Dataname, 
		     const blitz::Array<T, D>& Data)
{
  write_field(Dataname, Data, pred_data<T>());
}

template<int D> inline void
HdfFile::write_field(const std::string& Dataname, 
		     const blitz::Array<std::string, D>& Data)
{
  try {
  using namespace H5;
  create_group_if_needed(Dataname, *h);
  hsize_t dim[D];
  for(int i = 0; i < D; ++i)
    dim[i] = Data.extent(i);
  DataSpace ds(D, dim);
  StrType st(H5::PredType::C_S1, H5T_VARIABLE);
  DataSet d = h->createDataSet(Dataname, st, ds);
  blitz::Array<std::string, D> data2 = to_c_order_const(Data);
  blitz::Array<const char*, D> data3(data2.shape());
  typename blitz::Array<std::string, D>::const_iterator i1;
  typename blitz::Array<const char*, D>::iterator i2 = data3.begin();
  for(i1 = data2.begin(); i1 != data2.end(); ++i1, ++i2)
    *i2 = i1->c_str();
  d.write(data3.dataFirst(), st);
  write_attribute(Dataname + "/Type", "VarLenStr");
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While writing field " << Dataname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

template<int D> inline void
HdfFile::write_field(const std::string& Dataname, 
		     const blitz::Array<const char*, D>& Data)
{
  try {
  using namespace H5;

  // Figure out maximum length for char* since we 
  // are creating a fixed length string dataset
  int max_len = 0;
  typename blitz::Array<const char*, D>::const_iterator i1;
  for(i1 = Data.begin(); i1 != Data.end(); ++i1)
    max_len = std::max((int) strlen(*i1), max_len);
  // +1 for fairly odd convention of add a "\0" at the end
  H5::StrType dtype = H5::StrType(H5::PredType::C_S1, max_len + 1);

  create_group_if_needed(Dataname, *h);
  hsize_t dim[D];
  int flat_size = 1;
  for(int i = 0; i < D; ++i) {
    dim[i] = Data.extent(i);
    flat_size *= Data.extent(i);
  }

  DataSpace ds(D, dim);
  DataSet d = h->createDataSet(Dataname, dtype, ds);

  // Copy string data to a flattened version as expected
  // by HDF5 library
  blitz::Array<char, 2> data2(flat_size, max_len + 1);
  data2 = ' ';
  data2(blitz::Range::all(), max_len) = '\0';
  typename blitz::Array<char, 2>::iterator i2 = data2.begin();
  for(i1 = Data.begin(); i1 != Data.end(); ++i1) {
    strcpy(&(*i2), *i1);
    for(int j = 0; j < max_len + 1; j++)
      i2++;
  }

  d.write(data2.dataFirst(), dtype);
  write_type(Dataname, dtype);
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While writing field " << Dataname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

inline void HdfFile::write_attribute(const std::string& Aname, const char* Data)
{ write_attribute(Aname, std::string(Data)); }

//-----------------------------------------------------------------------
/// Write a given field. If you give a group (e.g., "Group1/Field1"),
/// then we automatically create the group if it doesn't already
/// exist.
///
/// This version lets you override the default type of the output data.
//-----------------------------------------------------------------------

template<class T, int D> inline void
HdfFile::write_field(const std::string& Dataname, 
		     const blitz::Array<T, D>& Data,
		     H5::DataType P)
{
  try {
  using namespace H5;
  create_group_if_needed(Dataname, *h);
  hsize_t dim[D];
  for(int i = 0; i < D; ++i)
    dim[i] = Data.extent(i);
  DataSpace ds(D, dim);
  DataSet d = h->createDataSet(Dataname, P, ds);
  blitz::Array<T, D> data2 = to_c_order_const(Data);
  d.write(data2.dataFirst(), pred_arr<T>());
  write_type(Dataname, P);
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While writing field " << Dataname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

//-----------------------------------------------------------------------
/// Write a given field. If you give a group (e.g., "Group1/Field1"),
/// then we automatically create the group if it doesn't already
/// exist.
///
/// This is a shortcut for writing an array that is exactly 1 element
/// long, this is the same as calling write_field with a
/// blitz::Array<T, 1> of size 1.
//-----------------------------------------------------------------------
template<class T> inline void HdfFile::write_field(const std::string& Dataname, 
					    const T& Data)
{
  write_field(Dataname, Data, pred_data<T>());
}

inline void HdfFile::write_field(const std::string& Dataname, 
				 const std::string& Data)
{
  try {
  using namespace H5;
  create_group_if_needed(Dataname, *h);
  hsize_t dim[1] = {1};
  DataSpace ds(1, dim);
  StrType st(PredType::C_S1, H5T_VARIABLE);
  DataSet d = h->createDataSet(Dataname, st, ds);
  const char* data2 = Data.c_str();
  d.write(&data2, st);
  write_attribute(Dataname + "/Type", "VarLenStr");
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While writing field " << Dataname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

inline void HdfFile::write_field(const std::string& Dataname, 
				 const char* Data)
{
  try {
  using namespace H5;
  create_group_if_needed(Dataname, *h);
  hsize_t dim[1] = {1};
  DataSpace ds(1, dim);
  // +1 is for the odd convention of including the trailing "\0" in
  // the data.
  H5::StrType dtype = H5::StrType(H5::PredType::C_S1, strlen(Data) + 1);
  DataSet d = h->createDataSet(Dataname, dtype, ds);
  d.write(Data, dtype);
  write_type(Dataname, dtype);
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While writing field " << Dataname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

//-----------------------------------------------------------------------
/// Write a given field. If you give a group (e.g., "Group1/Field1"),
/// then we automatically create the group if it doesn't already
/// exist.
///
/// This is a shortcut for writing an array that is exactly 1 element
/// long, this is the same as calling write_field with a
/// blitz::Array<T, 1> of size 1.
//-----------------------------------------------------------------------
template<class T> inline void HdfFile::write_field(const std::string& Dataname, 
					    const T& Data,
					    H5::DataType P)
{
  try {
  using namespace H5;
  create_group_if_needed(Dataname, *h);
  hsize_t dim[1] = {1};
  DataSpace ds(1, dim);
  DataSet d = h->createDataSet(Dataname, P, ds);
  d.write(&Data, pred_arr<T>());
  write_type(Dataname, P);
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While writing field " << Dataname
       << " for the file '" << fname
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}

}
#endif
