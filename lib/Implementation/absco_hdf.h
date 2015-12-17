#ifndef ABSCO_HDF_H
#define ABSCO_HDF_H
#include "absco.h"
#include "hdf_file.h"
#include "spectral_bound.h"
#include <blitz/array.h>
#include <vector>


namespace FullPhysics {
/****************************************************************//**
  This class is used to read the absco tables. This reads the HDF
  version of the files.

  Note that performance reasons we cache the data as we read it. The
  default cache is about 50 MB, which is a bit large but not too
  large. This can be adjusted if needed, either up for better
  performance or down for less memory.
*******************************************************************/
class AbscoHdf: public Absco {
public:
  AbscoHdf(const std::string& Fname, double Table_scale = 1.0, 
	   int Cache_nline = 5000);
  AbscoHdf(const std::string& Fname, 
	   const SpectralBound& Spectral_bound,
	   const std::vector<double>& Table_scale,
	   int Cache_nline = 5000);
  void load_file(const std::string& Fname);
  void load_file(const std::string& Fname, double Table_scale,
		 int Cache_nline = 5000);
  void load_file(const std::string& Fname, 
		 const SpectralBound& Spectral_bound,
		 const std::vector<double>& Table_scale,
		 int Cache_nline = 5000);
  virtual ~AbscoHdf() {}
  virtual std::string broadener_name() const { return bname; }
  virtual blitz::Array<double, 1> broadener_vmr_grid() const
  { return bvmr; }
  virtual double table_scale(double wn) const;
  virtual blitz::Array<double, 1> pressure_grid() const {return pgrid;}
  virtual blitz::Array<double, 2> temperature_grid() const {return tgrid;}
  virtual bool have_data(double wn) const;
  virtual bool is_float() const { return is_float_;}
  virtual std::string file_name() const { return hfile->file_name(); } 
  virtual void print(std::ostream& Os) const;
protected:
  virtual blitz::Array<double, 3> read_double(double wn) const;
  virtual blitz::Array<float, 3> read_float(double wn) const;
private:
  bool is_float_;
  int cache_nline;
  blitz::Array<double, 1> bvmr;
  boost::shared_ptr<HdfFile> hfile;
  mutable int cache_double_lbound;
  mutable int cache_double_ubound;
  mutable int cache_float_lbound;
  mutable int cache_float_ubound;
  template<class T> void bound_set(int lbound, int sz) const;
  mutable blitz::Array<double, 4> read_cache_double;
  mutable blitz::Array<float, 4> read_cache_float;
  template<class T> blitz::Array<T, 4>& read_cache() const;
  template<class T> void swap(int i) const;
  int wn_index(double Wn_in) const;
  std::string field_name;
  std::string bname;
  blitz::Array<double, 1> pgrid;
  blitz::Array<double, 2> tgrid;
  blitz::Array<double, 1> wn;
  SpectralBound sb;
  std::vector<double> table_scale_;
  // This just points to front and back of wn. We store this here so
  // don't need to keep calculating this.
  double *wnfront, *wnback;
};

template<> inline blitz::Array<double, 4>& AbscoHdf::read_cache<double>() const
{ return read_cache_double; }

template<> inline blitz::Array<float, 4>& AbscoHdf::read_cache<float>() const
{ return read_cache_float; }

template<> inline void AbscoHdf::bound_set<double>(int lbound, int sz) const
{
  cache_double_lbound = lbound;
  cache_double_ubound = lbound + sz;
}

template<> inline void AbscoHdf::bound_set<float>(int lbound, int sz) const
{
  cache_float_lbound = lbound;
  cache_float_ubound = lbound + sz;
}

}
#endif
