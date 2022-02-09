#include "ground_brdf.h"
#include "polynomial_eval.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

// Number of parameters in spars for the RT
#define NUM_BRDF_PARAMS 5

extern "C" {
    double black_sky_albedo_veg_f(const double* params, const double* sza);
    double black_sky_albedo_soil_f(const double* params, const double* sza);
    double exact_brdf_value_veg_f(const double* params, const double* sza, const double* vza, const double* azm);
    double exact_brdf_value_soil_f(const double* params, const double* sza, const double* vza, const double* azm);
}

#ifdef HAVE_LUA
#include "register_lua.h"

double black_sky_albedo_simple_veg(const blitz::Array<double, 1>& params, const double sza) {
    return black_sky_albedo_veg_f(params.dataFirst(), &sza);
}

double black_sky_albedo_simple_soil(const blitz::Array<double, 1>& params, const double sza) {
    return black_sky_albedo_soil_f(params.dataFirst(), &sza);
}

double exact_brdf_value_simple_veg(const blitz::Array<double, 1>& params, const double sza, const double vza, const double azm) {
    return exact_brdf_value_veg_f(params.dataFirst(), &sza, &vza, &azm);
}

double exact_brdf_value_simple_soil(const blitz::Array<double, 1>& params, const double sza, const double vza, const double azm) {
    return exact_brdf_value_soil_f(params.dataFirst(), &sza, &vza, &azm);
}

REGISTER_LUA_DERIVED_CLASS(GroundBrdfVeg, Ground)
.def(luabind::constructor<const blitz::Array<double, 2>&, const blitz::Array<bool, 2>&, const ArrayWithUnit<double, 1>&, const std::vector<std::string>&>())
.scope
[
    luabind::def("black_sky_albedo", &black_sky_albedo_simple_veg)
]
.scope
[
    luabind::def("kernel_value", &exact_brdf_value_simple_veg)
]
REGISTER_LUA_END()

REGISTER_LUA_DERIVED_CLASS(GroundBrdfSoil, Ground)
.def(luabind::constructor<const blitz::Array<double, 2>&, const blitz::Array<bool, 2>&, const ArrayWithUnit<double, 1>&, const std::vector<std::string>&>())
.scope
[
    luabind::def("black_sky_albedo", &black_sky_albedo_simple_soil)
]
.scope
[
    luabind::def("kernel_value", &exact_brdf_value_simple_soil)
]
REGISTER_LUA_END()
#endif

/****************************************************************//**
  Constructor that defines coefficients in a 2d array:
  Num_spectrometer * num_coeff
  Each row has the num_coeff BRDF parameters for a spectrometer.
  Coefficients are ordered:
  0: Rahman kernel factor
  1: Rahman hotspot parameter
  2: Rahman asymmetry factor
  3: Rahman anisotropy parameter
  4: Breon kernel factor
  5: Weight coefficent 1
     ...
  5+n-1: Weight coefficent n
 *******************************************************************/

GroundBrdf::GroundBrdf(const blitz::Array<double, 2>& Coeffs,
                       const blitz::Array<bool, 2>& Flag,
                       const ArrayWithUnit<double, 1>& Ref_points,
                       const std::vector<std::string>& Desc_band_names)
: reference_points(Ref_points), desc_band_names(Desc_band_names)
{
    if(Coeffs.cols() < NUM_BRDF_PARAMS + 1) {
        Exception err_msg;
        err_msg << "Number of parameters in Coeffs: " << Coeffs.cols() << " is < " << NUM_BRDF_PARAMS + 1 << " as expected";
        throw err_msg;
    }

    num_weight_params = Coeffs.cols() - NUM_BRDF_PARAMS;
    num_coeff = NUM_BRDF_PARAMS + num_weight_params;

    if(Coeffs.rows() != Flag.rows()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Coeffs: " << Coeffs.rows() << " does not match the number in Flag: " << Flag.rows();
        throw err_msg;
    }

    if(Coeffs.cols() != Flag.cols()) {
        Exception err_msg;
        err_msg << "Number of parameters in Coeffs: " << Coeffs.cols() << " does not match the number in Flag: " << Flag.cols();
        throw err_msg;
    }

    if(Coeffs.rows() != (int) Desc_band_names.size()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Coeffs: " << Coeffs.rows() << " does not match the number in Desc_band_names: " << Desc_band_names.size();
        throw err_msg;
    }

    if(Ref_points.rows() != (int) Desc_band_names.size()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Ref_points: " << Ref_points.rows() << " does not match the number in Desc_band_names: " << Desc_band_names.size();
        throw err_msg;
    }

    // Make local arrays to deal with const issues on call to init. The init routine copies the data
    Array<double, 2> coeffs(Coeffs);
    Array<bool, 2> flags(Flag);

    // Flatten arrays for state vector
    init(Array<double, 1>(coeffs.dataFirst(), TinyVector<int, 1>(Coeffs.rows() * Coeffs.cols()), blitz::neverDeleteData),
         Array<bool, 1>(flags.dataFirst(), TinyVector<int, 1>(Flag.rows() * Flag.cols()), blitz::neverDeleteData));
}

/// Protected constructor that matches the dimensionality of coeff and flag arrays
GroundBrdf::GroundBrdf(const blitz::Array<double, 1>& Spec_coeffs,
                       const blitz::Array<bool, 1>& Flag,
                       const ArrayWithUnit<double, 1>& Ref_points,
                       const std::vector<std::string>& Desc_band_names)
  : SubStateVectorArray<Ground>(Spec_coeffs, Flag),
    reference_points(Ref_points), desc_band_names(Desc_band_names)
{
    num_weight_params = Spec_coeffs.rows() / Ref_points.rows() - NUM_BRDF_PARAMS;
    num_coeff = NUM_BRDF_PARAMS + num_weight_params;
}

ArrayAd<double, 1> GroundBrdf::surface_parameter(const double wn, const int spec_index) const
{
    AutoDerivative<double> w = weight(wn, spec_index);
    ArrayAd<double, 1> spars;
    spars.resize(NUM_BRDF_PARAMS, coefficient().number_variable());
    spars(0) = w * rahman_factor(spec_index);
    spars(1) = hotspot_parameter(spec_index);
    spars(2) = asymmetry_parameter(spec_index);
    spars(3) = anisotropy_parameter(spec_index);
    spars(4) = w * breon_factor(spec_index);
    return spars;
}

const AutoDerivative<double> GroundBrdf::weight(const double wn, const int spec_index) const
{
    double ref_wn = reference_points(spec_index).convert_wave(units::inv_cm).value;
    ArrayAd<double, 1> weight_params(weight_parameters(spec_index));
    Poly1d weight_poly(weight_params, false);
    AutoDerivative<double> wn_ad(wn); // Make sure we use the AutoDerivative interface to Poly1d
    return weight_poly(wn_ad - ref_wn);
}

//----

const AutoDerivative<double> GroundBrdf::rahman_factor(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(num_coeff * spec_index + RAHMAN_KERNEL_FACTOR_INDEX);
}

const AutoDerivative<double> GroundBrdf::hotspot_parameter(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(num_coeff * spec_index + RAHMAN_OVERALL_AMPLITUDE_INDEX);
}

const AutoDerivative<double> GroundBrdf::asymmetry_parameter(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(num_coeff * spec_index + RAHMAN_ASYMMETRY_FACTOR_INDEX);
}

const AutoDerivative<double> GroundBrdf::anisotropy_parameter(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(num_coeff * spec_index + RAHMAN_GEOMETRIC_FACTOR_INDEX);
}

const AutoDerivative<double> GroundBrdf::breon_factor(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(num_coeff * spec_index + BREON_KERNEL_FACTOR_INDEX);
}

const AutoDerivative<double> GroundBrdf::weight_intercept(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(num_coeff * spec_index + BRDF_WEIGHT_INTERCEPT_INDEX);
}

const AutoDerivative<double> GroundBrdf::weight_slope(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(num_coeff * spec_index + BRDF_WEIGHT_SLOPE_INDEX);
}

AutoDerivative<double> GroundBrdf::weight_coeff(const int spec_index, const int weight_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(num_coeff * spec_index + (NUM_BRDF_PARAMS + weight_index));
}

ArrayAd<double, 1> GroundBrdf::weight_parameters(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    int offset = num_coeff * spec_index + NUM_BRDF_PARAMS;
    return coefficient()(Range(offset, offset + num_weight_params - 1));
}

//----

void GroundBrdf::rahman_factor(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(num_coeff * spec_index + RAHMAN_KERNEL_FACTOR_INDEX) = val;
}

void GroundBrdf::hotspot_parameter(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(num_coeff * spec_index + RAHMAN_OVERALL_AMPLITUDE_INDEX) = val;
}

void GroundBrdf::asymmetry_parameter(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(num_coeff * spec_index + RAHMAN_ASYMMETRY_FACTOR_INDEX) = val;
}

void GroundBrdf::anisotropy_parameter(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(num_coeff * spec_index + RAHMAN_GEOMETRIC_FACTOR_INDEX) = val;
}

void GroundBrdf::breon_factor(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(num_coeff * spec_index + BREON_KERNEL_FACTOR_INDEX) = val;
}

void GroundBrdf::weight_intercept(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(num_coeff * spec_index + BRDF_WEIGHT_INTERCEPT_INDEX) = val;
}

void GroundBrdf::weight_slope(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(num_coeff * spec_index + BRDF_WEIGHT_SLOPE_INDEX) = val;
}

void GroundBrdf::weight_coeff(const int spec_index, const AutoDerivative<double>& val, const int weight_index)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(num_coeff * spec_index + (NUM_BRDF_PARAMS + weight_index)) = val;
}

void GroundBrdf::weight_parameters(const int spec_index, const ArrayAd<double, 1>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    int offset = num_coeff * spec_index + NUM_BRDF_PARAMS;
    coeff(Range(offset, offset + num_weight_params - 1)) = val;
}

//----

const blitz::Array<double, 2> GroundBrdf::brdf_covariance(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    // Return empty array if covariance is empty due to not being retrieved
    if(product(statevector_covariance().shape()) == 0) {
        return Array<double, 2>(0, 0);
    }

    int ret_num_param = statevector_covariance().rows() / desc_band_names.size();
    int ret_cov_offset = ret_num_param * spec_index;

    blitz::Array<double, 2> cov(num_coeff, num_coeff);
    cov = 0.0;

    // Copy out the retrieved covariance values into a matrix for the spectrometer that includes
    // zeros for non retrieved elements
    int in_idx_a = ret_cov_offset;
    for (int out_idx_a = 0; out_idx_a < num_coeff; out_idx_a++) {
        if (used_flag_value()(out_idx_a)) {
            int in_idx_b = ret_cov_offset;
            for (int out_idx_b = 0; out_idx_b < num_coeff; out_idx_b++) {
                if (used_flag_value()(out_idx_b)) {
                    cov(out_idx_a, out_idx_b) = statevector_covariance()(in_idx_a, in_idx_b++);
                }
            }
            in_idx_a++;
        }
    }

    return cov;
}

// Helper function
blitz::Array<double, 1> GroundBrdf::black_sky_params(const int Spec_index)
{
    double ref_wn = reference_points(Spec_index).convert_wave(units::inv_cm).value;
    double w = weight(ref_wn, Spec_index).value();

    blitz::Array<double, 1> params(NUM_BRDF_PARAMS, blitz::ColumnMajorArray<1>());
    params(0) = w * rahman_factor(Spec_index).value();
    params(1) = hotspot_parameter(Spec_index).value();
    params(2) = asymmetry_parameter(Spec_index).value();
    params(3) = anisotropy_parameter(Spec_index).value();
    params(4) = w * breon_factor(Spec_index).value();
    return params;
}

// Helper function
blitz::Array<double, 1> GroundBrdf::kernel_value_params(const int Spec_index)
{
    blitz::Array<double, 1> params(NUM_BRDF_PARAMS, blitz::ColumnMajorArray<1>());
    params(0) = rahman_factor(Spec_index).value();
    params(1) = hotspot_parameter(Spec_index).value();
    params(2) = asymmetry_parameter(Spec_index).value();
    params(3) = anisotropy_parameter(Spec_index).value();
    params(4) = breon_factor(Spec_index).value();
    return params;
}

const double GroundBrdfVeg::black_sky_albedo(const int Spec_index, const double Sza)
{
    blitz::Array<double, 1> params = black_sky_params(Spec_index);
    return black_sky_albedo_veg_f(params.dataFirst(), &Sza);
}

const double GroundBrdfSoil::black_sky_albedo(const int Spec_index, const double Sza)
{
    blitz::Array<double, 1> params = black_sky_params(Spec_index);
    return black_sky_albedo_soil_f(params.dataFirst(), &Sza);
}

const double GroundBrdfVeg::kernel_value(const int Spec_index, const double Sza, const double Vza, const double Azm)
{
    blitz::Array<double, 1> params = kernel_value_params(Spec_index);
    return exact_brdf_value_veg_f(params.dataFirst(), &Sza, &Vza, &Azm);
}

const double GroundBrdfSoil::kernel_value(const int Spec_index, const double Sza, const double Vza, const double Azm)
{
    blitz::Array<double, 1> params = kernel_value_params(Spec_index);
    return exact_brdf_value_soil_f(params.dataFirst(), &Sza, &Vza, &Azm);
}

std::string GroundBrdf::state_vector_name_i(int i) const {
    int b_idx = int(i / num_coeff);
    int c_idx = i - num_coeff * b_idx;

    std::stringstream name;
    name << "Ground BRDF " << breon_type() << " " << desc_band_names[b_idx] << " ";
    if (c_idx == RAHMAN_KERNEL_FACTOR_INDEX)
        name << "Rahman Factor";
    else if (c_idx == RAHMAN_OVERALL_AMPLITUDE_INDEX)
        name << "Hotspot Parameter";
    else if (c_idx == RAHMAN_ASYMMETRY_FACTOR_INDEX)
        name << "Asymmetry Parameter";
    else if (c_idx == RAHMAN_GEOMETRIC_FACTOR_INDEX)
        name << "Anisotropy Parameter";
    else if (c_idx == BREON_KERNEL_FACTOR_INDEX)
        name << "Breon Factor";
    else if (c_idx >= NUM_BRDF_PARAMS && c_idx < num_coeff)
        name << "Weight coefficent " << c_idx - NUM_BRDF_PARAMS + 1;
    else
        name << "Unknown Index " << i;

    return name.str();
}

void GroundBrdf::print(std::ostream& Os) const
{
    Os << "GroundBrdf:\n";
    for(int spec_idx = 0; spec_idx < number_spectrometer(); spec_idx++) {
        Os << "    " << desc_band_names[spec_idx] << ":" << std::endl;
        OstreamPad opad(Os, "        ");
        opad << "Rahman Factor: " << rahman_factor(spec_idx).value() << std::endl
             << "Hotspot Parameter: " << hotspot_parameter(spec_idx).value() << std::endl
             << "Asymmetry Parameter: " << asymmetry_parameter(spec_idx).value() << std::endl
             << "Anisotropy Parameter: " << anisotropy_parameter(spec_idx).value() << std::endl
             << "Breon Factor: " << breon_factor(spec_idx).value() << std::endl;

             ArrayAd<double, 1> weight_params(weight_parameters(spec_idx));

             for(int i = 0; i < num_weight_params; i++)
                  opad << "Weight coefficent " << i << ": " << weight_params(i).value() << std::endl;

        opad.strict_sync();
    }
}
