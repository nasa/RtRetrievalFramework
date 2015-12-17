#include <model_measure_oco.h>


using namespace FullPhysics;
using namespace blitz;


ModelMeasureOCO::ModelMeasureOCO(const boost::shared_ptr<ForwardModel>& fm)
  : FM(fm),
    meas_units(fm->measured_radiance_all().spectral_range().units())
{}


void ModelMeasureOCO::model_eval()
{
  if(M.size() <= 0)
    radiance_from_fm(true);
}


void ModelMeasureOCO::jacobian_eval()
{
  if(K.size() <= 0)
    radiance_from_fm(false);
}


void ModelMeasureOCO::model_jacobian_eval()
{
  if((K.size() <= 0) or (M.size() <= 0))
    radiance_from_fm(false);
}

int ModelMeasureOCO::expected_parameter_size() const
{
  return FM->state_vector()->observer_claimed_size();
}


void ModelMeasureOCO::radiance_from_fm(bool Skip_jacobian)
{
  assert_parameter_set_correctly();
  FM->state_vector()->update_state(X);

  //  TEMPORARY
  //
  //  There should be no need for this error check block.
  //  The problem is that with our current model code 
  //  measurement may change.
  //
  Array<double, 1> temp_msrmnt(FM->measured_radiance_all().spectral_range().data());
  if(msrmnt.rows() != temp_msrmnt.rows())
    throw Exception("Measurement has changed during the retrieval. :( ");
  double msrmnt_L1_norm = sum(abs(msrmnt));
  if(msrmnt_L1_norm <= 0.0)
    throw Exception("Measurement has no signal (just 0). :( ");
  if( sum(abs(temp_msrmnt-msrmnt))/msrmnt_L1_norm > 0.0000001 )
    throw Exception("Measurement has changed during the retrieval. :( ");

  Spectrum rad_spec = FM->radiance_all(Skip_jacobian);
  SpectralRange rad_mod = rad_spec.spectral_range().convert(meas_units);
  M.reference(rad_mod.data_ad().value());
  assert_model_correct(M);
  M.makeUnique();
  if(!Skip_jacobian) {
    K.reference(rad_mod.data_ad().jacobian());
    assert_jacobian_correct(K);

    //  TEMPORARY
    //
    // Should go away after we end support for 
    // fixed pressure level grid.
    vanishing_params_update();

    K.makeUnique();
  }
}




//  TEMPORARY
//
// Should go away after we end support for 
// fixed pressure level grid.
void ModelMeasureOCO::vanishing_params_update()
{
  if(K.size() <= 0) return;

  //  Even with fixed-pressure-level the following 
  //  code block should not be necessary if the
  //  forward-model class is already setting to zero
  //  the columns of the Jacobian associated with
  //  the unused parameters (if any).
  //
  Array<bool, 1> used(FM->state_vector()->used_flag());
  if(used.rows() != K.cols())
    throw Exception("Columns of Jacobian and elements of used-flag not equal in numbers! :( ");
  for(int i=0; i<used.rows(); i++)
    if(!used(i)) K(Range::all(),i) = 0.0;
}


void ModelMeasureOCO::parameters(const blitz::Array<double, 1>& x)
{
  //  This check is needed because it is not obvious 
  //  whether or not the Forward Model does caching.
  //
  if(!parameters_different(x)) return;
  ModelMeasure::parameters(x);
  FM->state_vector()->update_state(x);
}
