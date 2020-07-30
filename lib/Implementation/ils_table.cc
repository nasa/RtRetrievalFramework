#include "ils_table.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(IlsTableLinear, IlsFunction)
.def(luabind::constructor<const HdfFile&,
                          const int,
                          const std::string&, const std::string&,
                          const std::string&>())
.def(luabind::constructor<const HdfFile&,
                          const int,
                          const std::string&, const std::string&>())
.def(luabind::constructor<const blitz::Array<double, 1>&,
                          const blitz::Array<double, 2>&,
                          const blitz::Array<double, 2>&,
                          const blitz::Array<double, 1>&,
                          const blitz::Array<bool, 1>&,
                          const std::string&, const std::string&,
                          bool>())
.def("band_name", &IlsTableLinear::band_name)
REGISTER_LUA_END()

REGISTER_LUA_DERIVED_CLASS(IlsTableLog, IlsFunction)
.def(luabind::constructor<const HdfFile&,
                          const int,
                          const std::string&, const std::string&,
                          const std::string&>())
.def(luabind::constructor<const HdfFile&,
                          const int,
                          const std::string&, const std::string&>())
.def(luabind::constructor<const blitz::Array<double, 1>&,
                          const blitz::Array<double, 2>&,
                          const blitz::Array<double, 2>&,
                          const blitz::Array<double, 1>&,
                          const blitz::Array<bool, 1>&,
                          const std::string&, const std::string&,
                          bool>())
.def("band_name", &IlsTableLog::band_name)
REGISTER_LUA_END()

#endif

//-----------------------------------------------------------------------
/// Constructor, which reads the ILS information from a HDF file. We
/// determine the HDF group to read as Hdf_group + "/ILS_" +
/// (spec_index + 1), e.g., "Instrument/ILS/ILS_1". We then read the
/// fields delta_lambda, wavenumber, response, band_name, hdf_band_name,
/// and function_type
//-----------------------------------------------------------------------

template<class Lint>
IlsTable<Lint>::IlsTable(const HdfFile& Hdf_static_input,
                         int Spec_index,
                         const std::string& Band_name,
                         const std::string& Hdf_band_name,
                         const std::string& Hdf_group)
: SubStateVectorArray<IlsFunction>(1.0, false),
  band_name_(Band_name), hdf_band_name_(Hdf_band_name),
  from_hdf_file(true),
  hdf_file_name(Hdf_static_input.file_name())
{
  hdf_group = Hdf_group + "/ILS_" + 
    boost::lexical_cast<std::string>(Spec_index + 1);
  Array<double, 1> wavenumber
    (Hdf_static_input.read_field<double, 1>(hdf_group + "/wavenumber"));
  Array<double, 2> delta_lambda
    (Hdf_static_input.read_field<double, 2>(hdf_group + "/delta_lambda"));
  Array<double, 2> response
    (Hdf_static_input.read_field<double, 2>(hdf_group + "/response"));
  std::string ftype = 
    Hdf_static_input.read_field<std::string>(hdf_group + "/function_type");
  if(ftype == "table")
    interpolate_wavenumber = false;
  else if(ftype == "interpol")
    interpolate_wavenumber = true;
  else {
    Exception e;
    e << "Unrecognized function_type given in the file " 
      << Hdf_static_input.file_name() << " HDF group " << hdf_group
      << ". The value was '" << ftype << "'. This should be one of "
      << "'table' or 'interpol'";
    throw e;
  }
  create_delta_lambda_to_response(wavenumber, delta_lambda, response);
}

template<class Lint>
void IlsTable<Lint>::create_delta_lambda_to_response
(const blitz::Array<double, 1>& Wavenumber, 
 const blitz::Array<double, 2>& Delta_lambda, 
 const blitz::Array<double, 2>& Response)
{
  wavenumber_.reference(Wavenumber);
  delta_lambda_.reference(Delta_lambda);
  response_.reference(Response);

  delta_lambda_to_response.clear();
  for(int i = 0; i < Wavenumber.rows(); ++i) {
    // Create AutoDerivative version of these arrays, for use by the
    // interpolator. 
    Array<ad, 1> 
      delta_lambda(arrad(Delta_lambda(i, Range::all())).to_array()),
      response(arrad(Response(i, Range::all())).to_array());
    delta_lambda_to_response[Wavenumber(i)] = 
      Lint(delta_lambda.begin(), delta_lambda.end(), response.begin(),
           LinearInterpolate<ad, ad>::OUT_OF_RANGE_EXTRAPOLATE);
  }
}

// See base class for description.
template<class Lint>
std::string IlsTable<Lint>::state_vector_name_i(int i) const
{
  std::string res = "ILS " + band_name_ + " Scale";
  return res;
}

// See base class for description
template<class Lint>
void IlsTable<Lint>::ils
(const AutoDerivative<double>& wn_center,
 const blitz::Array<double, 1>& wn,
 ArrayAd<double, 1>& res, bool jac_optimization) const
{
  // Note that this function is a bottleneck, so we have written this
  // for speed at the expense of some clarity.
  res.resize(wn.rows(), wn_center.number_variable());
  AutoDerivative<double> f;
  f.gradient().resize(wn_center.number_variable());
  AutoDerivative<double> t;
  t.gradient().resize(wn_center.number_variable());
  AutoDerivativeRef<double> tref(t.value(), t.gradient());

  // w will be (wn(i) - wn_center)/coeff(0). We do this calculation in 2 steps
  // to speed up the calculation.
  AutoDerivative<double> w;
  w.gradient().resize(wn_center.number_variable());
  AutoDerivative<double> scale;
  // See IlsFunction for description of this optimization.
  // We use d[wn_center,coeff]
  if(jac_optimization)
    scale = AutoDerivative<double>(coeff(0).value(), 1, 2);
  else
    scale = coeff(0);
  blitz::Array<double, 1> g2 = (-wn_center/scale).gradient();
  blitz::Array<double, 1> g1;
  if(!scale.is_constant())
    g1.reference((1/scale).gradient());
  w.gradient() = g2;

  // This returns the value for the first wn not less than wn_center.
  it inter = delta_lambda_to_response.lower_bound(wn_center.value());

  // By the convention, we want the value for the last wn less than
  // wn_center.
  if(inter != delta_lambda_to_response.begin())
    --inter;

  for(int i = 0; i < res.rows(); ++i) {
    w.value() = (wn(i) - wn_center.value())/scale.value();
    if(!scale.is_constant())
      w.gradient() = wn(i) * g1 + g2;
    inter->second.interpolate(w, res(i));
  }

  // Do an interpolation in the wavenumber direction, if requested.
  if(interpolate_wavenumber) {
    double wn1 = inter->first;
    ++inter;
    if(inter != delta_lambda_to_response.end()) {
      double wn2 = inter->first;
      f = (wn_center - wn1) / (wn2 - wn1);
      for(int i = 0; i < res.rows(); ++i) {
        w.value() = wn(i) - wn_center.value();
	w.value() = (wn(i) - wn_center.value())/scale.value();
	if(!scale.is_constant())
	  w.gradient() = wn(i) * g1 + g2;
        inter->second.interpolate(w, tref);
        double v = res(i).value();
        res(i).value_ref() = (1 - f.value()) * v + f.value() * t.value();
        res(i).gradient_ref() *= (1 - f.value());
        res(i).gradient_ref() += f.gradient() * (t.value() -  v) +
          f.value() * t.gradient();
      }
    }
  }
}

template<class Lint>
boost::shared_ptr<IlsFunction> IlsTable<Lint>::clone() const
{
  return boost::shared_ptr<IlsFunction>
    (new IlsTable(wavenumber_, delta_lambda_, response_, coeff.value(),
                  used_flag, band_name_, hdf_band_name_, interpolate_wavenumber));
}

template<class Lint>
void IlsTable<Lint>::print(std::ostream& Os) const 
{
  Os << "IlsTable for band " << band_name() << "\n"
     << "  Interpolate wavenumber: " << (interpolate_wavenumber ? "true" : "false")
     << "\n";
  if(from_hdf_file)
    Os << "  Hdf file name:          " << hdf_file_name << "\n"
       << "  Hdf group:              " << hdf_group;
  Os << "  Scale:                  (";
  for(int i = 0; i < coeff.rows() - 1; ++i)
    Os << coeff.value()(i) << ", ";
  Os << coeff.value()(coeff.rows() - 1) << ")\n"
     << "  Retrieve:               (";
  for(int i = 0; i < used_flag.rows() - 1; ++i)
    Os << (used_flag(i) ? "true" : "false") << ", ";
  Os << (used_flag(used_flag.rows() - 1) ? "true" : "false") << ")";
}
