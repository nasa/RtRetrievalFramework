#include "ground_lambertian.h"
#include "polynomial_eval.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(GroundLambertian, Ground)
.def(luabind::constructor<const blitz::Array<double, 2>&,
                          const blitz::Array<bool, 2>&, 
                          const ArrayWithUnit<double, 1>&,
                          const std::vector<std::string>&>())
REGISTER_LUA_END()
#endif

GroundLambertian::GroundLambertian(const blitz::Array<double, 2>& Spec_coeffs,
                                   const blitz::Array<bool, 2>& Flag, 
                                   const ArrayWithUnit<double, 1>& Ref_points,
                                   const std::vector<std::string>& Desc_band_names) 
: reference_points(Ref_points), desc_band_names(Desc_band_names)
{
    if(Spec_coeffs.rows() != Flag.rows()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Spec_coeffs: " << Spec_coeffs.rows() << " does not match the number in Flag: " << Flag.rows();
        throw err_msg;
    }

    if(Spec_coeffs.cols() != Flag.cols()) {
        Exception err_msg;
        err_msg << "Number of parameters in Spec_coeffs: " << Spec_coeffs.cols() << " does not match the number in Flag: " << Flag.cols();
        throw err_msg;
    }

    if(Spec_coeffs.rows() != Ref_points.rows()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Spec_coeffs: " << Spec_coeffs.rows() << " does not match the number in Ref_points: " << Ref_points.rows();
        throw err_msg;
    }

    if(Spec_coeffs.rows() != (int) Desc_band_names.size()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Spec_coeffs: " << Spec_coeffs.rows() << " does not match the number in Desc_band_names: " << Desc_band_names.size();
        throw err_msg;
    }

    // Make local arrays to deal with const issues on call to init. The init routine copies the data
    Array<double, 2> spec_coeffs(Spec_coeffs);
    Array<bool, 2> flags(Flag);

    // Flatten arrays for state vector
    init(Array<double, 1>(spec_coeffs.dataFirst(), TinyVector<int, 1>(Spec_coeffs.rows() * Spec_coeffs.cols()), blitz::neverDeleteData),
         Array<bool, 1>(flags.dataFirst(), TinyVector<int, 1>(Flag.rows() * Flag.cols()), blitz::neverDeleteData));
}

// Protected constructor that matches the dimensionality of coeff and flag arrays
GroundLambertian::GroundLambertian(const blitz::Array<double, 1>& Spec_coeffs,
                                   const blitz::Array<bool, 1>& Flag, 
                                   const ArrayWithUnit<double, 1>& Ref_points,
                                   const std::vector<std::string>& Desc_band_names)
  : SubStateVectorArray<Ground>(Spec_coeffs, Flag), 
    reference_points(Ref_points), desc_band_names(Desc_band_names) 
{
}

ArrayAd<double, 1> GroundLambertian::surface_parameter(const double wn, const int spec_index) const
{
    AutoDerivative<double> wn_albedo = albedo(DoubleWithUnit(wn, units::inv_cm), spec_index); 
    ArrayAd<double, 1> spars(1, wn_albedo.number_variable());
    spars(0) = wn_albedo;
    return spars; 
}

const AutoDerivative<double> GroundLambertian::albedo(const DoubleWithUnit Wave_point, const int spec_index) const
{
    // Evaluate in terms of wavenumber
    AutoDerivative<double> wn(Wave_point.convert_wave(units::inv_cm).value);
    double ref_wn = reference_points(spec_index).convert_wave(units::inv_cm).value;
    // Calculation of albedo from state vector value
    Poly1d surface_poly(albedo_coefficients(spec_index), false);
    return surface_poly(wn - ref_wn);
}

const ArrayAd<double, 1> GroundLambertian::albedo_coefficients(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    int offset = number_params() * spec_index;
    return coefficient()(Range(offset, offset + number_params() - 1));
}

const blitz::Array<double, 2> GroundLambertian::albedo_covariance(const int spec_index) const
{ 
    range_check(spec_index, 0, number_spectrometer());

    // Return empty array if covariance is empty due to not being retrieved
    if(product(statevector_covariance().shape()) == 0) {
        return Array<double, 2>(0, 0);
    }

    int num_params = number_params();
    int offset = num_params * spec_index;

    blitz::Array<double, 2> cov(num_params, num_params);
    for (int param_a = 0; param_a < num_params; param_a++) {
        for (int param_b = 0; param_b < num_params; param_b++) {
            cov(param_a, param_b) = statevector_covariance()(offset + param_a, offset + param_b);
        }
    }
    return cov;
}

boost::shared_ptr<Ground> GroundLambertian::clone() const {
    return boost::shared_ptr<Ground>(new GroundLambertian(coefficient().value(), used_flag_value(), reference_points, desc_band_names));
}

std::string GroundLambertian::state_vector_name_i(int i) const {
    int b_idx = int(i / number_params());
    int c_idx = i - number_params() * b_idx;
    return "Ground Lambertian " + desc_band_names[b_idx] + " Albedo Parm " + boost::lexical_cast<std::string>(c_idx + 1);
}

void GroundLambertian::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "GroundLambertian:" << std::endl;
    for(int b_idx = 0; b_idx < number_spectrometer(); b_idx++) {
        opad << "Band: " << desc_band_names[b_idx] << std::endl
             << "Coefficient: " << albedo_coefficients(b_idx).value() << std::endl
             << "Flag: ";
        for(int c_idx = 0; c_idx < number_params(); c_idx++) {
            opad << used_flag_value()(number_params() * b_idx + c_idx);
            if(c_idx < number_params()-1) {
                opad << ", ";
            }
        }
        opad << std::endl << "Reference Point: " << reference_points(b_idx) << std::endl;
    }
    opad.strict_sync();
}
