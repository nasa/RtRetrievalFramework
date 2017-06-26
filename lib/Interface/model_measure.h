#ifndef MODEL_MEASURE_H
#define MODEL_MEASURE_H
#include <model_state.h>


namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for models and measurements.
///
/// This interface class is designed for the purpose of fitting a 
/// parametrized mathematical model (a vector function) to the
/// measured data.  The class hierarchy rooted at this base class
/// implements three components that are needed for the fitting
/// process:
///   - The mathematical model (vector) and its Jacobian
///   - The measurement (vector) and its error covariance matrix
///   - The statistical analysis method for model fitting
///
/// The class hierarchy rooted at ModelMeasure is designed 
/// assuming only one of the following two statistical analysis
/// methods:
///   - Maximum Likelihood 
///   - Maximum A Posteriori
///
/// Assuming only the above two statistical analysis methods, this
/// class is designed to be the common interface needed to implement
/// the two statistical methods mentioned above.  Therefore, in the
/// light of some other statistical analysis method, that I don't
/// know about it, it may or may not be better to redesign this class
/// hierarchy in some other way.
///
/// The parameter handling capability of this class is inherited
/// from ProblemState class (through direct or indirect inheritance).
/// However, most real-life mathematical models and their Jacobians
/// are computationally expensive; therefore, this class inherits
/// ModelState from ProblemState class hierarchy to inherit the
/// capability of storing and possibly reusing the last evaluations
/// of the model and its Jacobian.
///
/// The measurement error covariance matrix is not necessarily a
/// diagonal matrix; however, in the implementation of this class
/// it is assumed to be a diagonal matrix.
//-----------------------------------------------------------------------

class ModelMeasure : public ModelState {

public:


//-----------------------------------------------------------------------
/// \brief Constructor
///
/// \param[in] measurement
///            The measurement to which the model is to be fitted
///
/// \param[in] measurement_error_cov
///            The measurement error covariance matrix assumed to
///            be diagonal
//-----------------------------------------------------------------------

  ModelMeasure(const blitz::Array<double, 1>& measurement, 
               const blitz::Array<double, 1>& measurement_error_cov);


//-----------------------------------------------------------------------
/// \brief Default Constructor
//-----------------------------------------------------------------------

  ModelMeasure() {}


  virtual ~ModelMeasure() {}


//-----------------------------------------------------------------------
/// \brief Evaluates the model at the currently set parameter values
///
/// This method must be implemented by the classes derived from 
/// this class.
///
/// The parameters (the point in the parameter space) must have
/// already been set before calling this method.  The parameters are
/// already set if one of the following methods is already called
/// successfully:
///   - parameters() (see ProblemState class)
///   - model_x() (see CostFunc class)
///   - jacobian_x()
///   - model_jacobian_x()
/// 
/// If the parameters are already set, then this method evaluate the
/// model at the currently set parameter values (point in the 
/// parameter space).
//-----------------------------------------------------------------------

  virtual void model_eval() = 0;


//-----------------------------------------------------------------------
/// \brief Evaluates and returns the model at the currently set
///        parameter values
///
/// All the comments on model_eval() method also apply to this method,
/// and in addition this method returns the evaluated model.
///
/// The size of the model vector can be obtained in advance
/// by calling measurement_size().
///
/// \return Evaluated model at the current point 
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> model()
  { model_eval(); return M.copy(); }


//-----------------------------------------------------------------------
/// \brief The model function with parameters
///
/// This method also evaluates the model; however, it sets the model
/// at the input new point and then evaluates the model.
///
/// The size of the model vector can be obtained in advance
/// by calling measurement_size().
///
/// \param[in] x
///            New set of parameters
///
/// \return Evaluated model at the input new point 
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> model_x(const blitz::Array<double, 1>& x)
  { parameters(x); return model(); }


//-----------------------------------------------------------------------
/// \brief Evaluates the Jacobian of the model at the currently set
///        parameter values
///
/// This method must be implemented by the classes derived from 
/// this class.
///
/// The parameters (the point in the parameter space) must have
/// already been set before calling this method.  The parameters are
/// already set if one of the following methods is already called
/// successfully:
///   - parameters() (see ProblemState class)
///   - model_x() (see CostFunc class)
///   - jacobian_x()
///   - model_jacobian_x()
/// 
/// If the parameters are already set, then this method evaluate the
/// Jacobian of the model at the currently set parameter values
/// (point in the parameter space).
//-----------------------------------------------------------------------

  virtual void jacobian_eval() = 0;


//-----------------------------------------------------------------------
/// \brief Evaluates and returns the Jacobian of the model at the
///        currently set parameter values
///
/// All the comments on jacobian_eval() method also apply to this
/// method, and in addition this method returns the evaluated Jacobian.
///
/// The sizes of the Jacobian matrix can be obtained in advance:
///   - For the number of its rows call measurement_size().
///   - For the number of its columns call expected_parameter_size()
///     (see ProblemState class).
///
/// \return Evaluated model Jacobian at the current point 
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> jacobian()
  { jacobian_eval(); return K.copy(); }


//-----------------------------------------------------------------------
/// \brief The model Jacobian with parameters
///
/// This method also evaluates the Jacobian; however, it sets the model
/// at the input new point and then evaluates its Jacobian.
///
/// The sizes of the Jacobian matrix can be obtained in advance:
///   - For the number of its rows call measurement_size().
///   - For the number of its columns call expected_parameter_size()
///     (see ProblemState class).
///
/// \param[in] x
///            New set of parameters
///
/// \return Evaluated model Jacobian at the input new point 
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> jacobian_x(const blitz::Array<double, 1>& x)
  { parameters(x); return jacobian(); }


//-----------------------------------------------------------------------
/// \brief Evaluates the model and its Jacobian at the currently set
///        parameter values
///
/// The parameters (the point in the parameter space) must have
/// already been set before calling this method.  The parameters are
/// already set if one of the following methods is already called
/// successfully:
///   - parameters() (see ProblemState class)
///   - model_x() (see CostFunc class)
///   - jacobian_x()
///   - model_jacobian_x()
/// 
/// If the parameters are already set, then this method evaluate the
/// model and its Jacobian at the currently set parameter values
/// (point in the parameter space).
//-----------------------------------------------------------------------

  virtual void model_jacobian_eval()
  { jacobian_eval(); model_eval(); }


//-----------------------------------------------------------------------
/// \brief Evaluates model and its Jacobian together
///
/// All the comments on model_jacobian_eval() method also apply to
/// this method, and in addition this method passes to the caller the
/// evaluated model and its Jacobian.
///
/// \param[out] m 
///             The model
///
/// \param[out] k
///             The Jacobian of the model
//-----------------------------------------------------------------------

  virtual void model_jacobian(
     blitz::Array<double, 1>& m, blitz::Array<double, 2>& k)
  { model_jacobian_eval(); m.reference(M.copy()); k.reference(K.copy()); }


//-----------------------------------------------------------------------
/// \brief Model and its Jacobian with parameters
///
/// This method passes to the caller the evaluated model and its
/// Jacobian after setting the problem at the input new point.
///
/// \param[in] x
///            New set of parameters
/// 
/// \param[out] m 
///             The model
///
/// \param[out] k
///             The Jacobian of the model
//-----------------------------------------------------------------------

  virtual void model_jacobian_x(const blitz::Array<double, 1>& x,
     blitz::Array<double, 1>& m, blitz::Array<double, 2>& k)
  { parameters(x); model_jacobian(m,k); }


//-----------------------------------------------------------------------
/// \brief Returns the measured data, to which the model is fit.
///
/// \return The measurement vector.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> measurement() const
  { return msrmnt.copy(); }


//-----------------------------------------------------------------------
/// \brief Returns the measurement error covariance
///        (implemented as a diagonal matrix).
///
/// \return The measurement error covariance.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> measurement_error_cov() const
  { return Se.copy(); }


//-----------------------------------------------------------------------
/// \brief Returns the size of the measurement data vector
///
/// This method returns the size of the measured data vector.
/// The following must also be equal to the value returned by
/// measurement_size():
///   - The size of the model data
///   - The number of the rows of the model Jacobian
///   - The number of the diagonal elements of the error
///     covariance matrix
///
/// \return The size of the measurement data vector
//-----------------------------------------------------------------------

  virtual int measurement_size() const
  { return msrmnt.rows(); }


//-----------------------------------------------------------------------
/// \brief Conditions that must be satisfied when a derived class
///        computes the model
///
/// This method is just the implementation of some conditions that
/// must be satisfied after a derived class computes the model.
/// The derived class itself calls this method to check the
/// computed model.
///
/// \param[in] m 
///            The computed model
//-----------------------------------------------------------------------

  virtual void assert_model_correct(const blitz::Array<double, 1>& m) const;


//-----------------------------------------------------------------------
/// \brief Conditions that must be satisfied when a derived class
///        computes the Jacobian of the model
///
/// This method is just the implementation of some conditions that
/// must be satisfied after a derived class computes the  Jacobian
/// of the model.  The derived class itself calls this method to 
/// check the computed Jacobian.
///
/// \param[in] k
///            The computed model
//-----------------------------------------------------------------------

  virtual void assert_jacobian_correct(const blitz::Array<double, 2>& k) const;


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "ModelMeasure"; }




//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// The following public methods are for convenience.
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
/// \brief Returns model and measurement difference
///        (model - measurement)
///
/// This method is for convenience.  It returns the difference of the
/// computed model and the measurement.  The difference is not called
/// residual on purpose.  The term residual will be used in the 
/// context of the Non-Linear (or Linear) Least Squares optimization.
///
/// Let the following be the computed model and the measured data
/// respectively
///   - M
///   - S
///
/// Then this method returns
/// \f[
///     M - S
/// \f]
///
/// \return Model and measurement difference (model - measurement)
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> model_measure_diff()
  { return blitz::Array<double, 1>(model() - measurement()); }


//-----------------------------------------------------------------------
/// \brief Returns model and measurement difference weighted by the
///        inverse of the Cholesky decomposition of the error
///        covariance matrix
///
/// This method is for convenience.  It returns the difference of the
/// computed model and measurement, model_measure_diff(), weighted by
/// the Cholesky decomposition (roughly speaking the square root) of 
/// the error covariance matrix.  In the name of this method "uncert"
/// is the short for uncertainty, but probably using the word 
/// uncertainty or uncert in the name of this method is a bad Idea.
/// In the implementation of this class the error covariance matrix
/// is implemented as a diagonal matrix.  However, if the matrix 
/// were not diagonal, still the purpose of this method would be to
/// weight the model and the measurement difference with the Cholesky
/// decomposition of the entire error covariance matrix.
///
/// Let the following be the computed model, the measured data, and
/// the measurement error covariance matrix respectively:
///   - M
///   - S
///   - Se
///
/// Then the Cholesky decomposition of the error covariance matrix is 
/// \f[
///     S_e = C_e C_e^T
/// \f]
/// and this method returns
/// \f[
///     C_e^{-1}(M-S)
/// \f]
///
/// \return Model and measurement difference weighted by the inverse
///         of the Cholesky decomposition of the error covariance
///         matrix
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> uncert_weighted_model_measure_diff()
  { return blitz::Array<double, 1>(model_measure_diff()/Se_chol); }


//-----------------------------------------------------------------------
/// \brief Returns the model Jacobian weighted by the inverse of the
///        Cholesky decomposition of the error covariance matrix
///
/// This method is for convenience.  It returns the model Jacobian
/// weighted by the Cholesky decomposition (roughly speaking the 
/// square root) of the error covariance matrix.
///
/// Let the following be the computed model Jacobian and
/// the measurement error covariance matrix respectively:
///   - K
///   - Se
///
/// Then the Cholesky decomposition of the error covariance matrix is 
/// \f[
///     S_e = C_e C_e^T
/// \f]
/// and this method returns
/// \f[
///     C_e^{-1}K
/// \f]
///
/// The method uncert_weighted_model_measure_diff() is another 
/// function of the parameters, and the method 
/// uncert_weighted_jacobian() is the Jacobian of 
/// uncert_weighted_model_measure_diff().
///
/// \return Model Jacobian weighted by the inverse of the Cholesky
///         decomposition of the error covariance matrix
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> uncert_weighted_jacobian();


//-----------------------------------------------------------------------
/// \brief Returns the inner product of the matrix returned by
///        the method uncert_weighted_jacobian() by itself.
///
/// Let the following be the computed model Jacobian and
/// the measurement error covariance matrix respectively:
///   - K
///   - Se
///
/// Then this method returns
/// \f[
///     K^TS_e^{-1}K
/// \f]
///
/// \return The inner product of the matrix returned by
///         the method uncert_weighted_jacobian() by itself.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> uncert_weighted_jac_inner_product();


protected:

  blitz::Array<double, 1> msrmnt;
  blitz::Array<double, 1> Se;

  // For convenience
  blitz::Array<double, 1> Se_chol;

};
}
#endif
