#include "register_lua.h"
#include <blitz/array.h>
using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Register various Blitz array types into Lua
//-----------------------------------------------------------------------

// double 1d
double blitz_double_array_1d_read(const Array<double, 1>& V, int i)
{ range_check(i, 0, V.rows()); return V(i); }
blitz::Array<double, 1> blitz_double_array_1d_slice(const Array<double, 1>& V, Range i)
{ return V(i); }
void blitz_double_array_1d_set_i(Array<double, 1>& V, int i, double v)
{ range_check(i, 0, V.rows()); V(i) = v; }
void blitz_double_array_1d_set_r_val(Array<double, 1>& V, Range r, double v)
{ V(r) = v; }
void blitz_double_array_1d_set_r_arr(Array<double, 1>& V, Range r, Array<double, 1>& v)
{ V(r) = v; }

std::string blitz_double_array_1d_tostring(const Array<double, 1>& V)
{
  std::ostringstream os;
  os << "Array<double, 1>: " << V;
  return os.str();
}

// int 1d
int blitz_int_array_1d_read(const Array<int, 1>& V, int i)
{ range_check(i, 0, V.rows()); return V(i); }
void blitz_int_array_1d_set_i(Array<int, 1>& V, int i, int v)
{ range_check(i, 0, V.rows()); V(i) = v; }
void blitz_int_array_1d_set_r(Array<int, 1>& V, Range r, int v)
{ V(r) = v; }

std::string blitz_int_array_1d_tostring(const Array<int, 1>& V)
{
  std::ostringstream os;
  os << "Array<int, 1>: " << V;
  return os.str();
}

// bool 1d
bool blitz_bool_array_1d_read(const Array<bool, 1>& V, int i)
{ range_check(i, 0, V.rows()); return V(i); }
void blitz_bool_array_1d_set_i(Array<bool, 1>& V, int i, bool v)
{ range_check(i, 0, V.rows()); V(i) = v; }
void blitz_bool_array_1d_set_r(Array<bool, 1>& V, Range r, bool v)
{ V(r) = v; }

std::string blitz_bool_array_1d_tostring(const Array<bool, 1>& V)
{
  std::ostringstream os;
  os << "Array<bool, 1>: " << V;
  return os.str();
}

// double 2d
double blitz_double_array_2d_read(const Array<double, 2>& V, int i, int j)
{ range_check(i, 0, V.rows()); 
  range_check(j, 0, V.cols()); 
  return V(i, j); }
blitz::Array<double, 2> blitz_double_array_2d_slice_rr(const Array<double, 2>& V, Range i, Range j)
{ return V(i, j); }
blitz::Array<double, 1> blitz_double_array_2d_slice_ir(const Array<double, 2>& V, int i, Range j)
{ range_check(i, 0, V.rows()); return V(i, j); }
blitz::Array<double, 1> blitz_double_array_2d_slice_ri(const Array<double, 2>& V, Range i, int j)
{ range_check(j, 0, V.cols()); return V(i, j); }

void blitz_double_array_2d_set_ii(Array<double, 2>& V, int i, int j, double v)
{ range_check(i, 0, V.rows()); 
  range_check(j, 0, V.cols()); 
  V(i, j) = v; }
void blitz_double_array_2d_set_ri(Array<double, 2>& V, Range i, int j, double v)
{ range_check(j, 0, V.cols()); 
  V(i, j) = v; }
void blitz_double_array_2d_set_ir(Array<double, 2>& V, int i, Range j, double v)
{ range_check(i, 0, V.rows()); V(i, j) = v; }
void blitz_double_array_2d_set_rr(Array<double, 2>& V, Range i, Range j, double v)
{ V(i, j) = v; }
void blitz_double_array_2d_set_arr_ri(Array<double, 2>& V, Range i, int j, Array<double, 1> v)
{ range_check(j, 0, V.cols()); V(i, j) = v; }
void blitz_double_array_2d_set_arr_ir(Array<double, 2>& V, int i, Range j, Array<double, 1> v)
{ range_check(i, 0, V.rows()); V(i, j) = v; }
void blitz_double_array_2d_set_arr_rr(Array<double, 2>& V, Range i, Range j, Array<double, 2> v)
{ V(i, j) = v; }

std::string blitz_double_array_2d_tostring(const Array<double, 2>& V)
{
  std::ostringstream os;
  os << "Array<double, 2>: " << V;
  return os.str();
}

// int 2d
int blitz_int_array_2d_read(const Array<int, 2>& V, int i, int j)
{ range_check(i, 0, V.rows()); 
  range_check(j, 0, V.cols()); 
  return V(i, j); }
blitz::Array<int, 2> blitz_int_array_2d_slice_rr(const Array<int, 2>& V, Range i, Range j)
{ return V(i, j); }
blitz::Array<int, 1> blitz_int_array_2d_slice_ir(const Array<int, 2>& V, int i, Range j)
{ range_check(i, 0, V.rows()); 
  return V(i, j); }
blitz::Array<int, 1> blitz_int_array_2d_slice_ri(const Array<int, 2>& V, Range i, int j)
{ range_check(j, 0, V.cols()); 
  return V(i, j); }

void blitz_int_array_2d_set_ii(Array<int, 2>& V, int i, int j, int v)
{ range_check(i, 0, V.rows()); 
  range_check(j, 0, V.cols()); 
  V(i, j) = v; }
void blitz_int_array_2d_set_ri(Array<int, 2>& V, Range i, int j, int v)
{range_check(j, 0, V.cols()); 
  V(i, j) = v; }
void blitz_int_array_2d_set_ir(Array<int, 2>& V, int i, Range j, int v)
{ range_check(i, 0, V.rows()); V(i, j) = v; }
void blitz_int_array_2d_set_rr(Array<int, 2>& V, Range i, Range j, int v)
{ V(i, j) = v; }
void blitz_int_array_2d_set_arr_ri(Array<int, 2>& V, Range i, int j, Array<int, 1> v)
{ range_check(j, 0, V.cols()); V(i, j) = v; }
void blitz_int_array_2d_set_arr_ir(Array<int, 2>& V, int i, Range j, Array<int, 1> v)
{ range_check(i, 0, V.rows()); V(i, j) = v; }
void blitz_int_array_2d_set_arr_rr(Array<int, 2>& V, Range i, Range j, Array<int, 2> v)
{ V(i, j) = v; }

std::string blitz_int_array_2d_tostring(const Array<int, 2>& V)
{
  std::ostringstream os;
  os << "Array<int, 2>: " << V;
  return os.str();
}

// bool 2d
bool blitz_bool_array_2d_read(const Array<bool, 2>& V, int i, int j)
{ range_check(i, 0, V.rows()); 
  range_check(j, 0, V.cols()); 
  return V(i, j); }
blitz::Array<bool, 2> blitz_bool_array_2d_slice_rr(const Array<bool, 2>& V, Range i, Range j)
{ return V(i, j); }
blitz::Array<bool, 1> blitz_bool_array_2d_slice_ir(const Array<bool, 2>& V, int i, Range j)
{ range_check(i, 0, V.rows()); 
  return V(i, j); }
blitz::Array<bool, 1> blitz_bool_array_2d_slice_ri(const Array<bool, 2>& V, Range i, int j)
{ range_check(j, 0, V.cols()); 
  return V(i, j); }

void blitz_bool_array_2d_set_ii(Array<bool, 2>& V, int i, int j, bool v)
{ range_check(i, 0, V.rows()); 
  range_check(j, 0, V.cols()); 
  V(i, j) = v; }
void blitz_bool_array_2d_set_ri(Array<bool, 2>& V, Range i, int j, bool v)
{ range_check(j, 0, V.cols()); 
  V(i, j) = v; }
void blitz_bool_array_2d_set_ir(Array<bool, 2>& V, int i, Range j, bool v)
{ range_check(i, 0, V.rows()); 
  V(i, j) = v; }
void blitz_bool_array_2d_set_rr(Array<bool, 2>& V, Range i, Range j, bool v)
{ V(i, j) = v; }
void blitz_bool_array_2d_set_arr_ri(Array<bool, 2>& V, Range i, int j, Array<bool, 1> v)
{ range_check(j, 0, V.cols()); V(i, j) = v; }
void blitz_bool_array_2d_set_arr_ir(Array<bool, 2>& V, int i, Range j, Array<bool, 1> v)
{ range_check(i, 0, V.rows()); V(i, j) = v; }
void blitz_bool_array_2d_set_arr_rr(Array<bool, 2>& V, Range i, Range j, Array<bool, 2> v)
{ V(i, j) = v; }

std::string blitz_bool_array_2d_tostring(const Array<bool, 2>& V)
{
  std::ostringstream os;
  os << "Array<bool, 2>: " << V;
  return os.str();
}

// double 3d
double blitz_double_array_3d_read(const Array<double, 3>& V, int i, int j, int k)
{ return V(i, j, k); }

blitz::Array<double, 3> blitz_double_array_3d_slice_rrr(const Array<double, 3>& V, Range i, Range j, Range k)
{ return V(i, j, k); }
blitz::Array<double, 2> blitz_double_array_3d_slice_irr(const Array<double, 3>& V, int i, Range j, Range k)
{ return V(i, j, k); }
blitz::Array<double, 2> blitz_double_array_3d_slice_rir(const Array<double, 3>& V, Range i, int j, Range k)
{ return V(i, j, k); }
blitz::Array<double, 2> blitz_double_array_3d_slice_rri(const Array<double, 3>& V, Range i, Range j, int k)
{ return V(i, j, k); }

blitz::Array<double, 1> blitz_double_array_3d_slice_rii(const Array<double, 3>& V, Range i, int j, int k)
{ return V(i, j, k); }
blitz::Array<double, 1> blitz_double_array_3d_slice_iri(const Array<double, 3>& V, int i, Range j, int k)
{ return V(i, j, k); }
blitz::Array<double, 1> blitz_double_array_3d_slice_iir(const Array<double, 3>& V, int i, int j, Range k)
{ return V(i, j, k); }


void blitz_double_array_3d_set_rrr(Array<double, 3>& V, Range i, Range j, Range k, double v)
{ V(i, j, k) = v; }

void blitz_double_array_3d_set_irr(Array<double, 3>& V, int i, Range j, Range k, double v)
{ V(i, j, k) = v; }
void blitz_double_array_3d_set_rir(Array<double, 3>& V, Range i, int j, Range k, double v)
{ V(i, j, k) = v; }
void blitz_double_array_3d_set_rri(Array<double, 3>& V, Range i, Range j, int k, double v)
{ V(i, j, k) = v; }

void blitz_double_array_3d_set_iir(Array<double, 3>& V, int i, int j, Range k, double v)
{ V(i, j, k) = v; }
void blitz_double_array_3d_set_iri(Array<double, 3>& V, int i, Range j, int k, double v)
{ V(i, j, k) = v; }
void blitz_double_array_3d_set_rii(Array<double, 3>& V, Range i, int j, int k, double v)
{ V(i, j, k) = v; }

void blitz_double_array_3d_set_iii(Array<double, 3>& V, int i, int j, int k, double v)
{ V(i, j, k) = v; }

void blitz_double_array_3d_set_arr_rrr(Array<double, 3>& V, Range i, Range j, Range k, Array<double, 3> v)
{ V(i, j, k) = v; }

void blitz_double_array_3d_set_arr_irr(Array<double, 3>& V, int i, Range j, Range k, Array<double, 2> v)
{ V(i, j, k) = v; }
void blitz_double_array_3d_set_arr_rir(Array<double, 3>& V, Range i, int j, Range k, Array<double, 2> v)
{ V(i, j, k) = v; }
void blitz_double_array_3d_set_arr_rri(Array<double, 3>& V, Range i, Range j, int k, Array<double, 2> v)
{ V(i, j, k) = v; }

void blitz_double_array_3d_set_arr_iir(Array<double, 3>& V, int i, int j, Range k, Array<double, 1> v)
{ V(i, j, k) = v; }
void blitz_double_array_3d_set_arr_iri(Array<double, 3>& V, int i, Range j, int k, Array<double, 1> v)
{ V(i, j, k) = v; }
void blitz_double_array_3d_set_arr_rii(Array<double, 3>& V, Range i, int j, int k, Array<double, 1> v)
{ V(i, j, k) = v; }

std::string blitz_double_array_3d_tostring(const Array<double, 3>& V)
{
  std::ostringstream os;
  os << "Array<double, 3>: "<< V;
  return os.str();
}

// int 3d
int blitz_int_array_3d_read(const Array<int, 3>& V, int i, int j, int k)
{ return V(i, j, k); }

blitz::Array<int, 3> blitz_int_array_3d_slice_rrr(const Array<int, 3>& V, Range i, Range j, Range k)
{ return V(i, j, k); }
blitz::Array<int, 2> blitz_int_array_3d_slice_irr(const Array<int, 3>& V, int i, Range j, Range k)
{ return V(i, j, k); }
blitz::Array<int, 2> blitz_int_array_3d_slice_rir(const Array<int, 3>& V, Range i, int j, Range k)
{ return V(i, j, k); }
blitz::Array<int, 2> blitz_int_array_3d_slice_rri(const Array<int, 3>& V, Range i, Range j, int k)
{ return V(i, j, k); }

blitz::Array<int, 1> blitz_int_array_3d_slice_rii(const Array<int, 3>& V, Range i, int j, int k)
{ return V(i, j, k); }
blitz::Array<int, 1> blitz_int_array_3d_slice_iri(const Array<int, 3>& V, int i, Range j, int k)
{ return V(i, j, k); }
blitz::Array<int, 1> blitz_int_array_3d_slice_iir(const Array<int, 3>& V, int i, int j, Range k)
{ return V(i, j, k); }


void blitz_int_array_3d_set_rrr(Array<int, 3>& V, Range i, Range j, Range k, int v)
{ V(i, j, k) = v; }

void blitz_int_array_3d_set_irr(Array<int, 3>& V, int i, Range j, Range k, int v)
{ V(i, j, k) = v; }
void blitz_int_array_3d_set_rir(Array<int, 3>& V, Range i, int j, Range k, int v)
{ V(i, j, k) = v; }
void blitz_int_array_3d_set_rri(Array<int, 3>& V, Range i, Range j, int k, int v)
{ V(i, j, k) = v; }

void blitz_int_array_3d_set_iir(Array<int, 3>& V, int i, int j, Range k, int v)
{ V(i, j, k) = v; }
void blitz_int_array_3d_set_iri(Array<int, 3>& V, int i, Range j, int k, int v)
{ V(i, j, k) = v; }
void blitz_int_array_3d_set_rii(Array<int, 3>& V, Range i, int j, int k, int v)
{ V(i, j, k) = v; }

void blitz_int_array_3d_set_iii(Array<int, 3>& V, int i, int j, int k, int v)
{ V(i, j, k) = v; }

void blitz_int_array_3d_set_arr_rrr(Array<int, 3>& V, Range i, Range j, Range k, Array<int, 3> v)
{ V(i, j, k) = v; }

void blitz_int_array_3d_set_arr_irr(Array<int, 3>& V, int i, Range j, Range k, Array<int, 2> v)
{ V(i, j, k) = v; }
void blitz_int_array_3d_set_arr_rir(Array<int, 3>& V, Range i, int j, Range k, Array<int, 2> v)
{ V(i, j, k) = v; }
void blitz_int_array_3d_set_arr_rri(Array<int, 3>& V, Range i, Range j, int k, Array<int, 2> v)
{ V(i, j, k) = v; }

void blitz_int_array_3d_set_arr_iir(Array<int, 3>& V, int i, int j, Range k, Array<int, 1> v)
{ V(i, j, k) = v; }
void blitz_int_array_3d_set_arr_iri(Array<int, 3>& V, int i, Range j, int k, Array<int, 1> v)
{ V(i, j, k) = v; }
void blitz_int_array_3d_set_arr_rii(Array<int, 3>& V, Range i, int j, int k, Array<int, 1> v)
{ V(i, j, k) = v; }

std::string blitz_int_array_3d_tostring(const Array<int, 3>& V)
{
  std::ostringstream os;
  os << "Array<int, 3>: "<< V;
  return os.str();
}

// bool 3d
bool blitz_bool_array_3d_read(const Array<bool, 1>& V, int i, int j, int k)
{ return V(i, j, k); }

void blitz_bool_array_3d_set_iii(Array<bool, 3>& V, int i, int j, int k, bool v)
{ V(i, j, k) = v; }
void blitz_bool_array_3d_set_iir(Array<bool, 3>& V, int i, int j, Range k, bool v)
{ V(i, j, k) = v; }
void blitz_bool_array_3d_set_iri(Array<bool, 3>& V, int i, Range j, int k, bool v)
{ V(i, j, k) = v; }
void blitz_bool_array_3d_set_irr(Array<bool, 3>& V, int i, Range j, Range k, bool v)
{ V(i, j, k) = v; }
void blitz_bool_array_3d_set_rii(Array<bool, 3>& V, Range i, int j, int k, bool v)
{ V(i, j, k) = v; }
void blitz_bool_array_3d_set_rir(Array<bool, 3>& V, Range i, int j, Range k, bool v)
{ V(i, j, k) = v; }
void blitz_bool_array_3d_set_rri(Array<bool, 3>& V, Range i, Range j, int k, bool v)
{ V(i, j, k) = v; }
void blitz_bool_array_3d_set_rrr(Array<bool, 3>& V, Range i, Range j, Range k, bool v)
{ V(i, j, k) = v; }

std::string blitz_bool_array_3d_tostring(const Array<bool, 3>& V)
{
  std::ostringstream os;
  os << "Array<bool, 3>: " << V;
  return os.str();
}

// double 4d
int blitz_double_array_4d_extent_4(const Array<double, 4>& V) {
  return V.extent(fourthDim);
}

double blitz_double_array_4d_read(const Array<double, 4>& V, int i, int j, int k, int l)
{ return V(i, j, k, l); }

blitz::Array<double, 4> blitz_double_array_4d_slice_rrrr(const Array<double, 4>& V, Range i, Range j, Range k, Range l)
{ return V(i, j, k, l); }

blitz::Array<double, 3> blitz_double_array_4d_slice_irrr(const Array<double, 4>& V, int i, Range j, Range k, Range l)
{ return V(i, j, k, l); }
blitz::Array<double, 3> blitz_double_array_4d_slice_rirr(const Array<double, 4>& V, Range i, int j, Range k, Range l)
{ return V(i, j, k, l); }
blitz::Array<double, 3> blitz_double_array_4d_slice_rrir(const Array<double, 4>& V, Range i, Range j, int k, Range l)
{ return V(i, j, k, l); }
blitz::Array<double, 3> blitz_double_array_4d_slice_rrri(const Array<double, 4>& V, Range i, Range j, Range k, int l)
{ return V(i, j, k, l); }

blitz::Array<double, 2> blitz_double_array_4d_slice_iirr(const Array<double, 4>& V, int i, int j, Range k, Range l)
{ return V(i, j, k, l); }
blitz::Array<double, 2> blitz_double_array_4d_slice_riir(const Array<double, 4>& V, Range i, int j, int k, Range l)
{ return V(i, j, k, l); }
blitz::Array<double, 2> blitz_double_array_4d_slice_rrii(const Array<double, 4>& V, Range i, Range j, int k, int l)
{ return V(i, j, k, l); }
blitz::Array<double, 2> blitz_double_array_4d_slice_irir(const Array<double, 4>& V, int i, Range j, int k, Range l)
{ return V(i, j, k, l); }
blitz::Array<double, 2> blitz_double_array_4d_slice_irri(const Array<double, 4>& V, int i, Range j, Range k, int l)
{ return V(i, j, k, l); }
blitz::Array<double, 2> blitz_double_array_4d_slice_riri(const Array<double, 4>& V, Range i, int j, Range k, int l)
{ return V(i, j, k, l); }

blitz::Array<double, 1> blitz_double_array_4d_slice_iiir(const Array<double, 4>& V, int i, int j, int k, Range l)
{ return V(i, j, k, l); }
blitz::Array<double, 1> blitz_double_array_4d_slice_riii(const Array<double, 4>& V, Range i, int j, int k, int l)
{ return V(i, j, k, l); }

void blitz_double_array_4d_set_rrrr(Array<double, 4>& V, Range i, Range j, Range k, Range l, double v)
{ V(i, j, k, l) = v; }

void blitz_double_array_4d_set_irrr(Array<double, 4>& V, int i, Range j, Range k, Range l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_rirr(Array<double, 4>& V, Range i, int j, Range k, Range l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_rrir(Array<double, 4>& V, Range i, Range j, int k, Range l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_rrri(Array<double, 4>& V, Range i, Range j, Range k, int l, double v)
{ V(i, j, k, l) = v; }

void blitz_double_array_4d_set_iirr(Array<double, 4>& V, int i, int j, Range k, Range l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_riir(Array<double, 4>& V, Range i, int j, int k, Range l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_rrii(Array<double, 4>& V, Range i, Range j, int k, int l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_irir(Array<double, 4>& V, int i, Range j, int k, Range l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_irri(Array<double, 4>& V, int i, Range j, Range k, int l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_riri(Array<double, 4>& V, Range i, int j, Range k, int l, double v)
{ V(i, j, k, l) = v; }

void blitz_double_array_4d_set_riii(Array<double, 4>& V, Range i, int j, int k, int l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_irii(Array<double, 4>& V, int i, Range j, int k, int l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_iiri(Array<double, 4>& V, int i, int j, Range k, int l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_iiir(Array<double, 4>& V, int i, int j, int k, Range l, double v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_iiii(Array<double, 4>& V, int i, int j, int k, int l, double v)
{ V(i, j, k, l) = v; }

void blitz_double_array_4d_set_arr_irrr(Array<double, 4>& V, int i, Range j, Range k, Range l, Array<double, 3>& v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_arr_rirr(Array<double, 4>& V, Range i, int j, Range k, Range l, Array<double, 3>& v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_arr_rrir(Array<double, 4>& V, Range i, Range j, int k, Range l, Array<double, 3>& v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_arr_rrri(Array<double, 4>& V, Range i, Range j, Range k, int l, Array<double, 3>& v)
{ V(i, j, k, l) = v; }

void blitz_double_array_4d_set_arr_iirr(Array<double, 4>& V, int i, int j, Range k, Range l, Array<double, 2>& v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_arr_riir(Array<double, 4>& V, Range i, int j, int k, Range l, Array<double, 2>& v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_arr_rrii(Array<double, 4>& V, Range i, Range j, int k, int l, Array<double, 2>& v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_arr_irir(Array<double, 4>& V, int i, Range j, int k, Range l, Array<double, 2>& v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_arr_irri(Array<double, 4>& V, int i, Range j, Range k, int l, Array<double, 2>& v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_arr_riri(Array<double, 4>& V, Range i, int j, Range k, int l, Array<double, 2>& v)
{ V(i, j, k, l) = v; }

void blitz_double_array_4d_set_arr_riii(Array<double, 4>& V, Range i, int j, int k, int l, Array<double, 1>& v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_arr_irii(Array<double, 4>& V, int i, Range j, int k, int l, Array<double, 1>& v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_arr_iiri(Array<double, 4>& V, int i, int j, Range k, int l, Array<double, 1>& v)
{ V(i, j, k, l) = v; }
void blitz_double_array_4d_set_arr_iiir(Array<double, 4>& V, int i, int j, int k, Range l, Array<double, 1>& v)
{ V(i, j, k, l) = v; }

std::string blitz_double_array_4d_tostring(const Array<double, 4>& V)
{
  std::ostringstream os;
  os << "Array<double, 4>: "<< V;
  return os.str();
}

// Math operations

template <class Vt>
Vt blitz_oper_add_double(const Vt& Arr, const double v) {
  Vt res( Arr + v );
  return res;
}

template <class Vt>
Vt double_oper_add_blitz(const double v, const Vt& Arr) {
  Vt res( v + Arr );
  return res;
}

template <class Vt>
Vt blitz_oper_add_array(const Vt& Arr, const Vt v) {
  Vt res( Arr + v );
  return res;
}

template <class Vt>
Vt blitz_oper_sub_double(const Vt& Arr, const double v) {
  Vt res( Arr - v );
  return res;
}

template <class Vt>
Vt double_oper_sub_blitz(const double v, const Vt& Arr) {
  Vt res( v - Arr );
  return res;
}

template <class Vt>
Vt blitz_oper_sub_array(const Vt& Arr, const Vt v) {
  Vt res( Arr - v );
  return res;
}

template <class Vt>
Vt blitz_oper_mul_double(const Vt& Arr, const double v) {
  Vt res( Arr * v );
  return res;
}

template <class Vt>
Vt double_oper_mul_blitz(const double v, const Vt& Arr) {
  Vt res( v * Arr );
  return res;
}

template <class Vt>
Vt blitz_oper_mul_array(const Vt& Arr, const Vt v) {
  Vt res( Arr * v );
  return res;
}

template <class Vt>
Vt blitz_oper_div_double(const Vt& Arr, const double v) {
  Vt res( Arr / v );
  return res;
}

template <class Vt>
Vt double_oper_div_blitz(const double v, const Vt& Arr) {
  Vt res( v / Arr );
  return res;
}

template <class Vt>
Vt blitz_oper_div_array(const Vt& Arr, const Vt& v) {
  Vt res( Arr / v );
  return res;
}

template <class Vt>
double blitz_sum(const Vt& Arr) {
  return sum(Arr);
}

template <class Vt>
double blitz_product(const Vt& Arr) {
  return product(Arr);
}

template <class Vt>
double blitz_mean(const Vt& Arr) {
  return mean(Arr);
}

template <class Vt>
double blitz_min(const Vt& Arr) {
  return min(Arr);
}

template <class Vt>
double blitz_max(const Vt& Arr) {
  return max(Arr);
}

template <class Vt>
Vt blitz_log(const Vt& Arr) {
  Vt res( log(Arr) );
  return res;
}

template <class Vt>
Vt blitz_exp(const Vt& Arr) {
  Vt res( exp(Arr) );
  return res;
}

typedef Array<double, 1> ad1;
typedef Array<double, 2> ad2;
typedef Array<double, 3> ad3;
typedef Array<double, 4> ad4;
typedef Array<bool, 1> ab1;
typedef Array<bool, 2> ab2;
typedef Array<bool, 3> ab3;
typedef Array<int, 1> ai1;
typedef Array<int, 2> ai2;
typedef Array<int, 3> ai3;
REGISTER_LUA_CLASS_NAME(ad1, Blitz_double_array_1d)
.def(luabind::constructor<int>())
.def("rows", &Array<double, 1>::rows)
.def("__call", &blitz_double_array_1d_read)
.def("__call", &blitz_double_array_1d_slice)
.def("set", &blitz_double_array_1d_set_i)
.def("set", &blitz_double_array_1d_set_r_val)
.def("set", &blitz_double_array_1d_set_r_arr)
.def("__tostring", &blitz_double_array_1d_tostring)
.def("__add", &blitz_oper_add_double<ad1>)
.def("__add", &double_oper_add_blitz<ad1>)
.def("__add", &blitz_oper_add_array<ad1>)
.def("__sub", &blitz_oper_sub_double<ad1>)
.def("__sub", &double_oper_sub_blitz<ad1>)
.def("__sub", &blitz_oper_sub_array<ad1>)
.def("__mul", &blitz_oper_mul_double<ad1>)
.def("__mul", &double_oper_mul_blitz<ad1>)
.def("__mul", &blitz_oper_mul_array<ad1>)
.def("__div", &blitz_oper_div_double<ad1>)
.def("__div", &double_oper_div_blitz<ad1>)
.def("__div", &blitz_oper_div_array<ad1>)
.def("sum", &blitz_sum<ad1>)
.def("product", &blitz_product<ad1>)
.def("mean", &blitz_mean<ad1>)
.def("min", &blitz_min<ad1>)
.def("max", &blitz_max<ad1>)
.def("log", &blitz_log<ad1>)
.def("exp", &blitz_exp<ad1>)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(ai1, Blitz_int_array_1d)
.def(luabind::constructor<int>())
.def("rows", &Array<int, 1>::rows)
.def("__call", &blitz_int_array_1d_read)
.def("set", &blitz_int_array_1d_set_i)
.def("set", &blitz_int_array_1d_set_r)
.def("__tostring", &blitz_int_array_1d_tostring)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(ab1, Blitz_bool_array_1d)
.def(luabind::constructor<int>())
.def("rows", &Array<bool, 1>::rows)
.def("__call", &blitz_bool_array_1d_read)
.def("set", &blitz_bool_array_1d_set_i)
.def("set", &blitz_bool_array_1d_set_r)
.def("__tostring", &blitz_bool_array_1d_tostring)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(ad2, Blitz_double_array_2d)
.def(luabind::constructor<int,int>())
.def("rows", &Array<double, 2>::rows)
.def("cols", &Array<double, 2>::cols)
.def("__call", &blitz_double_array_2d_read)
.def("__call", &blitz_double_array_2d_slice_rr)
.def("__call", &blitz_double_array_2d_slice_ir)
.def("__call", &blitz_double_array_2d_slice_ri)
.def("set", &blitz_double_array_2d_set_ii)
.def("set", &blitz_double_array_2d_set_ri)
.def("set", &blitz_double_array_2d_set_ir)
.def("set", &blitz_double_array_2d_set_rr)
.def("set", &blitz_double_array_2d_set_arr_ri)
.def("set", &blitz_double_array_2d_set_arr_ir)
.def("set", &blitz_double_array_2d_set_arr_rr)
.def("__tostring", &blitz_double_array_2d_tostring)
.def("__add", &blitz_oper_add_double<ad2>)
.def("__add", &double_oper_add_blitz<ad2>)
.def("__add", &blitz_oper_add_array<ad2>)
.def("__sub", &blitz_oper_sub_double<ad2>)
.def("__sub", &double_oper_sub_blitz<ad2>)
.def("__sub", &blitz_oper_sub_array<ad2>)
.def("__mul", &blitz_oper_mul_double<ad2>)
.def("__mul", &double_oper_mul_blitz<ad2>)
.def("__mul", &blitz_oper_mul_array<ad2>)
.def("__div", &blitz_oper_div_double<ad2>)
.def("__div", &double_oper_div_blitz<ad2>)
.def("__div", &blitz_oper_div_array<ad2>)
.def("sum", &blitz_sum<ad2>)
.def("product", &blitz_product<ad2>)
.def("mean", &blitz_mean<ad2>)
.def("min", &blitz_min<ad2>)
.def("max", &blitz_max<ad2>)
.def("log", &blitz_log<ad2>)
.def("exp", &blitz_exp<ad2>)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(ai2, Blitz_int_array_2d)
.def(luabind::constructor<int,int>())
.def("rows", &Array<int, 2>::rows)
.def("cols", &Array<int, 2>::cols)
.def("__call", &blitz_int_array_2d_read)
.def("__call", &blitz_int_array_2d_slice_rr)
.def("__call", &blitz_int_array_2d_slice_ir)
.def("__call", &blitz_int_array_2d_slice_ri)
.def("set", &blitz_int_array_2d_set_ii)
.def("set", &blitz_int_array_2d_set_ri)
.def("set", &blitz_int_array_2d_set_ir)
.def("set", &blitz_int_array_2d_set_rr)
.def("set", &blitz_int_array_2d_set_arr_ri)
.def("set", &blitz_int_array_2d_set_arr_ir)
.def("set", &blitz_int_array_2d_set_arr_rr)
.def("__tostring", &blitz_int_array_2d_tostring)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(ab2, Blitz_bool_array_2d)
.def(luabind::constructor<int, int>())
.def("rows", &Array<bool, 2>::rows)
.def("cols", &Array<bool, 2>::cols)
.def("__call", &blitz_bool_array_2d_read)
.def("__call", &blitz_bool_array_2d_slice_rr)
.def("__call", &blitz_bool_array_2d_slice_ir)
.def("__call", &blitz_bool_array_2d_slice_ri)
.def("set", &blitz_bool_array_2d_set_ii)
.def("set", &blitz_bool_array_2d_set_ri)
.def("set", &blitz_bool_array_2d_set_ir)
.def("set", &blitz_bool_array_2d_set_rr)
.def("set", &blitz_bool_array_2d_set_arr_ri)
.def("set", &blitz_bool_array_2d_set_arr_ir)
.def("set", &blitz_bool_array_2d_set_arr_rr)
.def("__tostring", &blitz_bool_array_2d_tostring)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(ad3, Blitz_double_array_3d)
.def(luabind::constructor<int,int,int>())
.def("rows", &Array<double, 3>::rows)
.def("cols", &Array<double, 3>::cols)
.def("depth", &Array<double, 3>::depth)
.def("__call", &blitz_double_array_3d_read)
.def("__call", &blitz_double_array_3d_slice_rrr)
.def("__call", &blitz_double_array_3d_slice_irr)
.def("__call", &blitz_double_array_3d_slice_rir)
.def("__call", &blitz_double_array_3d_slice_rri)
.def("__call", &blitz_double_array_3d_slice_rii)
.def("__call", &blitz_double_array_3d_slice_iri)
.def("__call", &blitz_double_array_3d_slice_iir)
.def("set", &blitz_double_array_3d_set_iii)
.def("set", &blitz_double_array_3d_set_iir)
.def("set", &blitz_double_array_3d_set_iri)
.def("set", &blitz_double_array_3d_set_irr)
.def("set", &blitz_double_array_3d_set_rii)
.def("set", &blitz_double_array_3d_set_rir)
.def("set", &blitz_double_array_3d_set_rri)
.def("set", &blitz_double_array_3d_set_rrr)
.def("set", &blitz_double_array_3d_set_arr_iir)
.def("set", &blitz_double_array_3d_set_arr_iri)
.def("set", &blitz_double_array_3d_set_arr_irr)
.def("set", &blitz_double_array_3d_set_arr_rii)
.def("set", &blitz_double_array_3d_set_arr_rir)
.def("set", &blitz_double_array_3d_set_arr_rri)
.def("set", &blitz_double_array_3d_set_arr_rrr)
.def("__tostring", &blitz_double_array_3d_tostring)
.def("__add", &blitz_oper_add_double<ad3>)
.def("__add", &double_oper_add_blitz<ad3>)
.def("__add", &blitz_oper_add_array<ad3>)
.def("__sub", &blitz_oper_sub_double<ad3>)
.def("__sub", &double_oper_sub_blitz<ad3>)
.def("__sub", &blitz_oper_sub_array<ad3>)
.def("__mul", &double_oper_mul_blitz<ad3>)
.def("__mul", &double_oper_mul_blitz<ad3>)
.def("__mul", &blitz_oper_mul_array<ad3>)
.def("__div", &blitz_oper_div_double<ad3>)
.def("__div", &double_oper_div_blitz<ad3>)
.def("__div", &blitz_oper_div_array<ad3>)
.def("sum", &blitz_sum<ad3>)
.def("product", &blitz_product<ad3>)
.def("mean", &blitz_mean<ad3>)
.def("min", &blitz_min<ad3>)
.def("max", &blitz_max<ad3>)
.def("log", &blitz_log<ad3>)
.def("exp", &blitz_exp<ad3>)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(ai3, Blitz_int_array_3d)
.def(luabind::constructor<int,int,int>())
.def("rows", &Array<int, 3>::rows)
.def("cols", &Array<int, 3>::cols)
.def("depth", &Array<int, 3>::depth)
.def("__call", &blitz_int_array_3d_read)
.def("__call", &blitz_int_array_3d_slice_rrr)
.def("__call", &blitz_int_array_3d_slice_irr)
.def("__call", &blitz_int_array_3d_slice_rir)
.def("__call", &blitz_int_array_3d_slice_rri)
.def("__call", &blitz_int_array_3d_slice_rii)
.def("__call", &blitz_int_array_3d_slice_iri)
.def("__call", &blitz_int_array_3d_slice_iir)
.def("set", &blitz_int_array_3d_set_iii)
.def("set", &blitz_int_array_3d_set_iir)
.def("set", &blitz_int_array_3d_set_iri)
.def("set", &blitz_int_array_3d_set_irr)
.def("set", &blitz_int_array_3d_set_rii)
.def("set", &blitz_int_array_3d_set_rir)
.def("set", &blitz_int_array_3d_set_rri)
.def("set", &blitz_int_array_3d_set_rrr)
.def("set", &blitz_int_array_3d_set_arr_iir)
.def("set", &blitz_int_array_3d_set_arr_iri)
.def("set", &blitz_int_array_3d_set_arr_irr)
.def("set", &blitz_int_array_3d_set_arr_rii)
.def("set", &blitz_int_array_3d_set_arr_rir)
.def("set", &blitz_int_array_3d_set_arr_rri)
.def("set", &blitz_int_array_3d_set_arr_rrr)
.def("__tostring", &blitz_int_array_3d_tostring)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(ab3, Blitz_bool_array_3d)
.def(luabind::constructor<int>())
.def("rows", &Array<bool, 3>::rows)
.def("cols", &Array<bool, 3>::cols)
.def("depth", &Array<bool, 3>::depth)
.def("__call", &blitz_bool_array_3d_read)
.def("set", &blitz_bool_array_3d_set_iii)
.def("set", &blitz_bool_array_3d_set_iir)
.def("set", &blitz_bool_array_3d_set_iri)
.def("set", &blitz_bool_array_3d_set_irr)
.def("set", &blitz_bool_array_3d_set_rii)
.def("set", &blitz_bool_array_3d_set_rir)
.def("set", &blitz_bool_array_3d_set_rri)
.def("set", &blitz_bool_array_3d_set_rrr)
.def("__tostring", &blitz_bool_array_3d_tostring)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(ad4, Blitz_double_array_4d)
.def(luabind::constructor<int,int,int,int>())
.def("rows", &Array<double, 4>::rows)
.def("cols", &Array<double, 4>::cols)
.def("depth", &Array<double, 4>::depth)
.def("fourth_dim", &blitz_double_array_4d_extent_4)
.def("__call", &blitz_double_array_4d_read)
.def("__call", &blitz_double_array_4d_slice_rrrr)
.def("__call", &blitz_double_array_4d_slice_irrr)
.def("__call", &blitz_double_array_4d_slice_rirr)
.def("__call", &blitz_double_array_4d_slice_rrir)
.def("__call", &blitz_double_array_4d_slice_rrri)
.def("__call", &blitz_double_array_4d_slice_iirr)
.def("__call", &blitz_double_array_4d_slice_riir)
.def("__call", &blitz_double_array_4d_slice_rrii)
.def("__call", &blitz_double_array_4d_slice_irir)
.def("__call", &blitz_double_array_4d_slice_irri)
.def("__call", &blitz_double_array_4d_slice_riri)
.def("__call", &blitz_double_array_4d_slice_iiir)
.def("__call", &blitz_double_array_4d_slice_riii)
.def("set", &blitz_double_array_4d_set_rrrr)
.def("set", &blitz_double_array_4d_set_irrr)
.def("set", &blitz_double_array_4d_set_rirr)
.def("set", &blitz_double_array_4d_set_rrir)
.def("set", &blitz_double_array_4d_set_rrri)
.def("set", &blitz_double_array_4d_set_iirr)
.def("set", &blitz_double_array_4d_set_riir)
.def("set", &blitz_double_array_4d_set_rrii)
.def("set", &blitz_double_array_4d_set_irir)
.def("set", &blitz_double_array_4d_set_irri)
.def("set", &blitz_double_array_4d_set_riri)
.def("set", &blitz_double_array_4d_set_riii)
.def("set", &blitz_double_array_4d_set_irii)
.def("set", &blitz_double_array_4d_set_iiri)
.def("set", &blitz_double_array_4d_set_iiir)
.def("set", &blitz_double_array_4d_set_iiii)
.def("set", &blitz_double_array_4d_set_arr_irrr)
.def("set", &blitz_double_array_4d_set_arr_rirr)
.def("set", &blitz_double_array_4d_set_arr_rrir)
.def("set", &blitz_double_array_4d_set_arr_rrri)
.def("set", &blitz_double_array_4d_set_arr_iirr)
.def("set", &blitz_double_array_4d_set_arr_riir)
.def("set", &blitz_double_array_4d_set_arr_rrii)
.def("set", &blitz_double_array_4d_set_arr_irir)
.def("set", &blitz_double_array_4d_set_arr_irri)
.def("set", &blitz_double_array_4d_set_arr_riri)
.def("set", &blitz_double_array_4d_set_arr_riii)
.def("set", &blitz_double_array_4d_set_arr_irii)
.def("set", &blitz_double_array_4d_set_arr_iiri)
.def("set", &blitz_double_array_4d_set_arr_iiir)
.def("__tostring", &blitz_double_array_4d_tostring)
.def("__add", &blitz_oper_add_double<ad4>)
.def("__add", &double_oper_add_blitz<ad4>)
.def("__add", &blitz_oper_add_array<ad4>)
.def("__sub", &blitz_oper_sub_double<ad4>)
.def("__sub", &double_oper_sub_blitz<ad4>)
.def("__sub", &blitz_oper_sub_array<ad4>)
.def("__mul", &blitz_oper_mul_double<ad4>)
.def("__mul", &double_oper_mul_blitz<ad4>)
.def("__mul", &blitz_oper_mul_array<ad4>)
.def("__div", &blitz_oper_div_double<ad4>)
.def("__div", &double_oper_div_blitz<ad4>)
.def("__div", &blitz_oper_div_array<ad4>)
.def("sum", &blitz_sum<ad4>)
.def("product", &blitz_product<ad4>)
.def("mean", &blitz_mean<ad4>)
.def("min", &blitz_min<ad4>)
.def("max", &blitz_max<ad4>)
.def("log", &blitz_log<ad4>)
.def("exp", &blitz_exp<ad4>)
REGISTER_LUA_END()


std::string range_tostring(const blitz::Range& R)
{
  std::ostringstream os;
  os << R;
  return os.str();
}

REGISTER_LUA_CLASS_NAME(Range, Range)
.def(luabind::constructor<>())
.def(luabind::constructor<int>())
.def(luabind::constructor<int,int>())
.def(luabind::constructor<int,int,int>())
.def("__tostring", &range_tostring)
.scope
[
 luabind::def("all", &Range::all)
 ]
.enum_("Range")
[
 luabind::value("fromStart", (int) blitz::fromStart),
 luabind::value("toEnd", (int) blitz::toEnd)
]
REGISTER_LUA_END()
