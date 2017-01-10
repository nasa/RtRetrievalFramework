#include "hdf_file.h"

using namespace FullPhysics;
using namespace blitz;
using namespace H5;

#ifdef HAVE_LUA
#include "register_lua.h"
double hdf_file_apriori_double(const HdfFile& h, const std::string& gname)
{
  return h.read_field<double, 1>(gname + "/a_priori")(0);
}

std::string hdf_file_read_string(const HdfFile& h, const std::string& fname)
{
  return h.read_field<std::string>(fname);
}

std::vector<std::string> hdf_file_read_string_vector(const HdfFile& h, const std::string& fname)
{
  blitz::Array<std::string, 1> str_arr = h.read_field<std::string, 1>(fname);
  std::vector<std::string> result;
  for(int i = 0; i < str_arr.rows(); i++)
    result.push_back(str_arr(i));
  return result;
}

std::vector<std::string> hdf_file_read_string_vector_row(const HdfFile& h, const std::string& fname, int row)
{
  blitz::Array<std::string, 2> str_arr = h.read_field<std::string, 2>(fname);
  std::vector<std::string> result;
  for(int i = 0; i < str_arr.cols(); i++)
    result.push_back(str_arr(0,i));
  return result;
}

Array<int, 1> hdf_file_read_int_1d(const HdfFile& h, const std::string& fname)
{
  return h.read_field<int, 1>(fname);
}

Array<int, 2> hdf_file_read_int_2d(const HdfFile& h, const std::string& fname)
{
  return h.read_field<int, 2>(fname);
}

Array<int, 3> hdf_file_read_int_3d(const HdfFile& h, const std::string& fname)
{
  return h.read_field<int, 3>(fname);
}

double hdf_file_read_double_scalar(const HdfFile& h, const std::string& gname)
{
  return h.read_field<double>(gname);
}

Array<double, 1> hdf_file_read_double_1d(const HdfFile& h, 
					 const std::string& fname)
{
  return h.read_field<double, 1>(fname);
}

double hdf_file_read_double_1d_i(const HdfFile& h, 
				 const std::string& fname,
				 int I)
{
  blitz::TinyVector<int,1> start;
  blitz::TinyVector<int,1> sz;  
  start = I;
  sz = 1;
  return h.read_field<double, 1>(fname, start, sz)(0);
}

ArrayWithUnit<double, 1> hdf_file_read_double_with_unit_1d(const HdfFile& h, 
							   const std::string& fname)
{
  return h.read_field_with_unit<double, 1>(fname);
}


Array<double, 2> hdf_file_read_double_2d(const HdfFile& h, 
					 const std::string& fname)
{
  return h.read_field<double, 2>(fname);
}

double hdf_file_read_double_2d_i(const HdfFile& h, 
			       const std::string& fname,
			       int I, int J)
{
  blitz::TinyVector<int,2> start;
  blitz::TinyVector<int,2> sz;  
  start = I, J;
  sz = 1, 1;
  return h.read_field<double, 2>(fname, start, sz)(0, 0);
}

ArrayWithUnit<double, 2> hdf_file_read_double_with_unit_2d(const HdfFile& h, 
							   const std::string& fname)
{
  return h.read_field_with_unit<double, 2>(fname);
}


Array<double, 3> hdf_file_read_double_3d(const HdfFile& h, 
					 const std::string& fname)
{
  return h.read_field<double, 3>(fname);
}

double hdf_file_read_double_3d_i(const HdfFile& h, 
			       const std::string& fname,
			       int I, int J, int K)
{
  blitz::TinyVector<int,3> start;
  blitz::TinyVector<int,3> sz;  
  start = I, J, K;
  sz = 1, 1, 1;
  return h.read_field<double, 3>(fname, start, sz)(0, 0, 0);
}

ArrayWithUnit<double, 3> hdf_file_read_double_with_unit_3d(const HdfFile& h, 
					 const std::string& fname)
{
  return h.read_field_with_unit<double, 3>(fname);
}

Array<double, 4> hdf_file_read_double_4d(const HdfFile& h, 
					 const std::string& fname)
{
  return h.read_field<double, 4>(fname);
}

// We have a number of places where we read a 4d field and only want
// one sounding number. Have code for doing this.
Array<double, 3> hdf_file_read_double_4d_sounding(const HdfFile& h,
						  const std::string& fname,
						  int Sounding_num)
{
  TinyVector<int, 4> shp = h.read_shape<4>(fname);
  TinyVector<int,4> start;
  TinyVector<int,4> sz;  
  start = 0, Sounding_num, 0, 0;
  sz = shp(0), 1, shp(2), shp(3);
  return h.read_field<double, 4>(fname, start, sz)(Range::all(), 0, 
						   Range::all(), Range::all());
}


double hdf_file_read_double_4d_i(const HdfFile& h, 
			       const std::string& fname,
			       int I, int J, int K, int L)
{
  blitz::TinyVector<int,4> start;
  blitz::TinyVector<int,4> sz;  
  start = I, J, K, L;
  sz = 1, 1, 1, 1;
  return h.read_field<double, 4>(fname, start, sz)(0, 0, 0, 0);
}

Array<double, 5> hdf_file_read_double_5d(const HdfFile& h, 
					 const std::string& fname)
{
  return h.read_field<double, 5>(fname);
}

Array<double, 1> hdf_file_apriori(const HdfFile& h, 
				  const std::string& gname)
{
  return h.read_field<double, 1>(gname + "/a_priori");
}

ArrayWithUnit<double, 1> hdf_file_apriori_with_unit(const HdfFile& h, 
						    const std::string& gname)
{
  return h.read_field_with_unit<double, 1>(gname + "/a_priori");
}

Array<double, 1> hdf_file_apriori2(const HdfFile& h, 
				   const std::string& gname,
				   int row)
{
  return h.read_field<double, 2>(gname + "/a_priori")(row, Range::all());
}

Array<double, 1> hdf_file_apriori3(const HdfFile& h, 
				   const std::string& gname,
				   int row, int col)
{
  return h.read_field<double, 3>(gname + "/a_priori")(row, col, Range::all());
}

ArrayWithUnit<double, 1> hdf_file_apriori2_with_unit(const HdfFile& h, 
						     const std::string& gname,
						     int row)
{
  ArrayWithUnit<double, 1> res = h.read_field_with_unit<double, 2>(gname + "/a_priori")(row, Range::all());
  return res;
}

Array<double, 2> hdf_file_covariance(const HdfFile& h, 
				     const std::string& gname)
{
  return h.read_field<double, 2>(gname + "/covariance");
}

Array<double, 2> hdf_file_covariance2(const HdfFile& h, 
				      const std::string& gname,
				      int row)
{
  return h.read_field<double, 3>(gname + "/covariance")(row, Range::all(), 
							Range::all());
}

Array<double, 2> hdf_file_covariance3(const HdfFile& h, 
				      const std::string& gname,
				      int row, int col)
{
  return h.read_field<double, 4>(gname + "/covariance")(row, col, Range::all(), 
							Range::all());
}

REGISTER_LUA_CLASS(HdfFile)
.def(luabind::constructor<std::string>())
.def("read_string", &hdf_file_read_string)
.def("read_string_vector", &hdf_file_read_string_vector)
.def("read_string_vector", &hdf_file_read_string_vector_row)
.def("read_int_1d", &hdf_file_read_int_1d)
.def("read_int_2d", &hdf_file_read_int_2d)
.def("read_int_3d", &hdf_file_read_int_3d)
.def("read_double_scalar", &hdf_file_read_double_scalar)
.def("read_double_1d", &hdf_file_read_double_1d)
.def("read_double_1d", &hdf_file_read_double_1d_i)
.def("read_double_with_unit_1d", &hdf_file_read_double_with_unit_1d)
.def("read_double_2d", &hdf_file_read_double_2d)
.def("read_double_2d", &hdf_file_read_double_2d_i)
.def("read_double_with_unit_2d", &hdf_file_read_double_with_unit_2d)
.def("read_double_3d", &hdf_file_read_double_3d)
.def("read_double_3d", &hdf_file_read_double_3d_i)
.def("read_double_with_unit_3d", &hdf_file_read_double_with_unit_3d)
.def("read_double_4d", &hdf_file_read_double_4d)
.def("read_double_4d_sounding", &hdf_file_read_double_4d_sounding)
.def("read_double_5d", &hdf_file_read_double_5d)
.def("apriori", &hdf_file_apriori)
.def("apriori_with_unit", &hdf_file_apriori_with_unit)
.def("apriori", &hdf_file_apriori2)
.def("apriori_with_unit", &hdf_file_apriori2_with_unit)
.def("apriori", &hdf_file_apriori3)
.def("apriori_double", &hdf_file_apriori_double)
.def("covariance", &hdf_file_covariance)
.def("covariance", &hdf_file_covariance2)
.def("covariance", &hdf_file_covariance3)
.def("has_object", &HdfFile::has_object)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Open the given file with the given mode.
//-----------------------------------------------------------------------

HdfFile::HdfFile(const std::string& Fname, Mode M)
: fname(Fname), mode_(M)
{
  unsigned int flag;
  switch(M) {
  case READ:
    flag = H5F_ACC_RDONLY;
    break;
  case READ_WRITE:
    flag = H5F_ACC_RDWR;
    break;
  case CREATE:
    flag = H5F_ACC_TRUNC;
    break;
  default:
    throw Exception("This shouldn't happen");
  }
  H5::Exception::dontPrint();	// We'll report exceptions ourself,
				// don't have HDF library print
				// Warning messages.
  try {
    h.reset(new H5File(Fname, flag));
  } catch(const H5::Exception& e) {
    Exception en;
    en << "While trying to open file '" << Fname 
       << "' a HDF 5 Exception thrown:\n"
       << "  " << e.getDetailMsg();
    throw en;
  }
}
    
//-----------------------------------------------------------------------
/// Determine if Dataname contains a group, and if it does then create
/// it if we haven't already done so.
//-----------------------------------------------------------------------

void HdfFile::create_group_if_needed(const std::string& Dataname, 
				     H5::CommonFG& Parent)
{
  size_t i = Dataname.find_first_of('/');
  if(i ==std::string::npos)
    return;
  if(i ==0)			// If "/" is at beginning, then just
				// strip it off
    create_group_if_needed(Dataname.substr(i + 1), Parent);
  else {
    std::string gname = Dataname.substr(0, i);
    if(!is_present(gname, Parent)) {
      Group g = Parent.createGroup(gname);
      create_group_if_needed(Dataname.substr(i + 1), g);
    } else {
      Group g = Parent.openGroup(gname);
      create_group_if_needed(Dataname.substr(i + 1), g);
    }
  }
}

//-----------------------------------------------------------------------
/// Determine of object is a group or dataset.
//-----------------------------------------------------------------------

bool HdfFile::is_group(const std::string& Objname) const
{
  H5G_stat_t statbuf;
  h->getObjinfo(Objname, statbuf);
  return statbuf.type == H5G_GROUP;
}

//-----------------------------------------------------------------------
/// Determine if object is present.
//-----------------------------------------------------------------------

bool HdfFile::is_present(const std::string& Objname, 
			 const H5::CommonFG& Parent) const
{
  try {
    H5G_stat_t statbuf;
    Parent.getObjinfo(Objname, statbuf);
  } catch(const H5::Exception& E) {
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------
/// Determine if attribute is present.
//-----------------------------------------------------------------------

bool HdfFile::has_attribute(const std::string& Aname) const
{
  try {
    open_attribute(Aname);
  } catch(const H5::Exception& E) {
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------
/// Open an attribute attached to a Group. This should have a name
/// like "Dimensions/Aerosol/Size"
//-----------------------------------------------------------------------

H5::Attribute HdfFile::open_attribute(const std::string& Aname) const
{
  size_t i = Aname.find_last_of('/');
  if(i == std::string::npos) {
    FullPhysics::Exception e;
    e << "The attribute name '" << Aname << "' is not a valid name";
    throw e;
  }
  std::string gname = Aname.substr(0, i);
  std::string aname = Aname.substr(i + 1);
  if(is_group(gname)) {
    Group attr_group = h->openGroup(gname);
    return attr_group.openAttribute(aname);
  } else {
    DataSet attr_dataset = h->openDataSet(gname);
    return attr_dataset.openAttribute(aname);
  }
}

//-----------------------------------------------------------------------
/// Read the units for a dataset.
//-----------------------------------------------------------------------

Unit HdfFile::read_units(const std::string& Dataname) const
{
  Unit res;
  // Older OCO data does not have the "Units" attribute, but rather
  // has a older "Unit" attribute. Provide limited support for this
  // older format, so we can read the older data. This support can go
  // away in the future when we don't need this support any longer.
  if(has_attribute(Dataname + "/Unit") &&
     !has_attribute(Dataname + "/Units")) {
    std::string u = read_attribute<std::string>(Dataname + "/Unit");
    if(u == "Rotation_deg")
      res = units::deg;
    else if(u == "Distance_m")
      res = units::m;
    else if(u == "Velocity")
      res = units::m / units::s;
    else if(u == "Radiance")
      res = units::W / pow(units::cm, 2) / units::sr / units::inv_cm;
    else if(u == "Wavelength_microns")
      res = units::micron;
    else
      throw Exception("Unrecognized unit '" + u + " found for attribute "
		      + Dataname + "/Unit");
  } else
    // Otherwise, just process the Units string.
    res = Unit(read_attribute<std::string>(Dataname + "/Units"));
  return res;
}

//-----------------------------------------------------------------------
/// Create a new attribute
//-----------------------------------------------------------------------

H5::Attribute HdfFile::create_attribute(const std::string& Aname, 
		const H5::DataSpace& Ds, const H5::DataType& P)
{
  size_t i = Aname.find_last_of('/');
  if(i == std::string::npos) {
    FullPhysics::Exception e;
    e << "The attribute name '" << Aname << "' is not a valid name";
    throw e;
  }
  std::string gname = Aname.substr(0, i);
  std::string aname = Aname.substr(i + 1);
  if(!is_present(gname, *h))
    create_group_if_needed(Aname, *h);
  if(is_group(gname)) {
    Group attr_group = h->openGroup(gname);
    return attr_group.createAttribute(aname, P, Ds);
  } else {
    DataSet attr_dataset = h->openDataSet(gname);
    return attr_dataset.createAttribute(aname, P, Ds);
  }
}

//-----------------------------------------------------------------------
/// The SDOS products have a particular metadata convention where the
/// information about each of the Dimensions that appear in the file
/// are documented. This function takes this metadata for dimension
/// and puts it into the file.
//-----------------------------------------------------------------------

void HdfFile::dimension_metadata(const std::string& Name,
				 const std::string& Description,
				 int Size)
{
  write_attribute("/Dimensions/" + Name + "/Size", Size);
  write_attribute("/Dimensions/" + Name + "/Description", Description);
}

//-----------------------------------------------------------------------
/// The SDOS products have a particular metadata convention where the
/// information about each of the Shapes that appear in the file
/// are documented. This function takes this metadata for Shape
/// and puts it into the file.
//-----------------------------------------------------------------------

void HdfFile::shape_metadata(const std::string& Name, const std::string& Dim1)
{
  write_attribute("/Shapes/" + Name + "/Rank", 1);
  write_attribute("/Shapes/" + Name + "/Dimensions", Dim1);
}

//-----------------------------------------------------------------------
/// The SDOS products have a particular metadata convention where the
/// information about each of the Shapes that appear in the file
/// are documented. This function takes this metadata for Shape
/// and puts it into the file.
//-----------------------------------------------------------------------

void HdfFile::shape_metadata(const std::string& Name, const std::string& Dim1,
			     const std::string& Dim2)
{
  write_attribute("/Shapes/" + Name + "/Rank", 2);
  std::vector<std::string> dim;
  dim.push_back(Dim1);
  dim.push_back(Dim2);
  write_attribute("/Shapes/" + Name + "/Dimensions", dim);
}

//-----------------------------------------------------------------------
/// The SDOS products have a particular metadata convention where the
/// information about each of the Shapes that appear in the file
/// are documented. This function takes this metadata for Shape
/// and puts it into the file.
//-----------------------------------------------------------------------

void HdfFile::shape_metadata(const std::string& Name, const std::string& Dim1,
			     const std::string& Dim2,
			     const std::string& Dim3)
{
  write_attribute("/Shapes/" + Name + "/Rank", 3);
  std::vector<std::string> dim;
  dim.push_back(Dim1);
  dim.push_back(Dim2);
  dim.push_back(Dim3);
  write_attribute("/Shapes/" + Name + "/Dimensions", dim);
}

//-----------------------------------------------------------------------
/// The SDOS products have a particular metadata convention where the
/// information about each of the types appear as metadata. This is a
/// bit redundant since the HDF file already marks the type, put
/// presumably in some context reading the string metadata is
/// easier. In any case, this writes out the type for the recognized
/// types. We just silently skip the metadata if we don't recognize
/// the type.
//-----------------------------------------------------------------------

void HdfFile::write_type(const std::string& Dataname, const H5::DataType& P)
{
  if(P ==PredType::NATIVE_INT32)
    write_attribute(Dataname + "/Type", "Int32");
  else if(P ==PredType::NATIVE_FLOAT)
    write_attribute(Dataname + "/Type", "Float32");
  // we just ignore this if we don't recognize the type.
}

void HdfFile::print(std::ostream& Os) const
{
  Os << "HdfFile: \n" 
     << "  File name: " << file_name() << "\n";
}
