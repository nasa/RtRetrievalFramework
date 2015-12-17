#ifndef MAX_A_POSTERIORI_H
#define MAX_A_POSTERIORI_H
#include <model_measure.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for maximum a posteriori estimation.
///
/// This class is the base class for all classes that use 
/// maximum a posteriori estimation method to implement the problem
/// of estimating the parameters of a statistical model.
//-----------------------------------------------------------------------

class MaxAPosteriori : 
    virtual public ModelMeasure {

public:


//-----------------------------------------------------------------------
/// \brief Constructor
///
/// \param[in] a_priori_params
///            A priori knowledge on the parameters
///
/// \param[in] a_priori_cov
///            A priori covariance matrix
//-----------------------------------------------------------------------

  MaxAPosteriori(const blitz::Array<double, 1>& a_priori_params,
                 const blitz::Array<double, 2>& a_priori_cov);


  virtual ~MaxAPosteriori() {}


//-----------------------------------------------------------------------
/// \brief Returns the a priori values (knowledge) of the parameters.
///
/// \return The a priori values (knowledge) of the parameters
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> a_priori_params() const
  { return Xa.copy(); }


//-----------------------------------------------------------------------
/// \brief Returns the a priori covariance matrix
///
/// \return The a priori covariance matrix
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> a_priori_cov() const
  { return Sa.copy(); }


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "MaxAPosteriori"; }


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// The following public methods are for convenience.
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
/// \brief Returns the current parameters value and their a priori
///        value difference (current param - a priori param)
///
/// This method is for convenience.  It returns the difference of the
/// current parameters value and their a priori value. 
///
/// Let the following be the current parameters value (a vector) and
/// the parameters a prior value (another vector) respectively
///   - X
///   - Xa
///
/// Then this method returns
/// \f[
///     X - X_a
/// \f]
///
/// \return The current parameters values and their a-priori
///         value difference (current param - a priori param)
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> parameter_a_priori_diff() const
  { return blitz::Array<double, 1>(X-Xa); }


//-----------------------------------------------------------------------
/// \brief Returns the current parameters value and their a priori
///        value difference (current param - a priori param) weighted
///        by the inverse of the Cholesky decomposition of the a priori
///        covariance matrix
///
/// This method is for convenience.  It returns the current
/// parameters value and their a priori value difference weighted
/// by the inverse of the Cholesky decomposition of the a priori
/// covariance matrix. 
///
/// This method does not have a good name. The "cov_weighted" portion
/// of its name suggests that the difference is weighted by the 
/// a priori covariance matrix, but it is weighted by the inverse of
/// the Cholesky decomposition of the covariance matrix.
///
/// Let the following be the current parameters value, the parameters
/// a priori value, and the a priori covariance matrix respectively:
///   - X
///   - Xa
///   - Sa
///
/// Then the Cholesky decomposition of the a priori covariance matrix is 
/// \f[
///     S_a = C_a C_a^T
/// \f]
/// and this method returns
/// \f[
///     C_a^{-1}(X-X_a)
/// \f]
///
/// \return The current parameters value and their a priori
///         value difference (current param - a priori param) weighted
///         by the inverse of the Cholesky decomposition of the a priori
///         covariance matrix
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> cov_weighted_parameter_a_priori_diff() const;


//-----------------------------------------------------------------------
/// \brief Returns the inverse of the Cholesky decomposition of the
///        a priori covariance matrix
///
/// Let the following be the a priori covariance matrix:
///   - Sa
///
/// Then the Cholesky decomposition of the a priori covariance matrix is 
/// \f[
///     S_a = C_a C_a^T
/// \f]
/// and this method returns
/// \f[
///     C_a^{-1}
/// \f]
///
/// The method cov_weighted_parameter_a_priori_diff() is the
/// implementation of a function of the parameters and the
/// method a_priori_cov_chol_inv() is the Jacobian of
/// cov_weighted_parameter_a_priori_diff().
///
/// \return The inverse of the Cholesky decomposition of the a priori
///         covariance matrix
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> a_priori_cov_chol_inv() const
  { return Sa_chol_inv.copy(); }


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// The following are secondary convenient public methods.
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
/// \brief Returns the vector returned by model_measure_diff()
///        augmented at the bottom by the vector returned by
///        parameter_a_priori_diff()
///
/// This method is for convenience.  It returns an augmented vector.
/// The vectors returned by the following methods appear at the top 
/// and at the bottom of the augmented vector respectively:
///   - model_measure_diff()
///   - parameter_a_priori_diff()
///
/// Let the following be the computed model, the measured data,
/// the current parameters value, and the parameters a priori value
/// respectively:
///   - M
///   - S
///   - X
///   - Xa
///
/// Then this method returns
/// \f[
///     \left[ \begin{array}{c}
///              M - S 
///              ---- 
///              X - X_a
///            \end{array} \right]
/// \f]
///
/// \return The vector returned by model_measure_diff()
///         augmented at the bottom by the vector returned by
///         parameter_a_priori_diff()
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> model_measure_diff_aug();


//-----------------------------------------------------------------------
/// \brief Returns the vector returned by 
///        uncert_weighted_model_measure_diff()
///        augmented at the bottom by the vector returned by
///        cov_weighted_parameter_a_priori_diff()
///
/// This method is for convenience.  It returns an augmented vector.
/// The vectors returned by the following methods appear at the top 
/// and at the bottom of the augmented vector respectively:
///   - uncert_weighted_model_measure_diff()
///   - cov_weighted_parameter_a_priori_diff()
///
/// Assume the following:
///   - M (computed model)
///   - S (measurement data)
///   - Se (measurement error covariance matrix)
///   - X (current parameters value)
///   - Xa (parameters a priori value)
///   - Sa (a priori covariance matrix)
///
/// The Cholesky decomposition of the error covariance matrix is 
/// \f[
///     S_e = C_e C_e^T
/// \f]
/// and the Cholesky decomposition of the a priori covariance matrix is 
/// \f[
///     S_a = C_a C_a^T
/// \f]
///
/// Then this method returns
/// \f[
///     \left[ \begin{array}{c}
///              C_e^{-1}(M-S) 
///              ------ 
///              C_a^{-1}(X-X_a)
///            \end{array} \right]
/// \f]
///
/// \return The vector returned by 
///         uncert_weighted_model_measure_diff()
///         augmented at the bottom by the vector returned by
///         cov_weighted_parameter_a_priori_diff()
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> weighted_model_measure_diff_aug();


//-----------------------------------------------------------------------
/// \brief Returns the matrix returned by 
///        uncert_weighted_jacobian()
///        augmented at the bottom by the matrix returned by
///        a_priori_cov_chol_inv()
///
/// This method is for convenience.  It returns an augmented matrix.
/// The matrices returned by the following methods appear at the top 
/// and at the bottom of the augmented matrix respectively:
///   - uncert_weighted_jacobian()
///   - a_priori_cov_chol_inv()
///
/// Assume the following:
///   - K (computed model Jacobian)
///   - Se (measurement error covariance matrix)
///   - Sa (a priori covariance matrix)
///
/// The Cholesky decomposition of the error covariance matrix is 
/// \f[
///     S_e = C_e C_e^T
/// \f]
/// and the Cholesky decomposition of the a priori covariance matrix is 
/// \f[
///     S_a = C_a C_a^T
/// \f]
///
/// Then this method returns
/// \f[
///     \left[ \begin{array}{c}
///              C_e^{-1}K 
///              --- 
///              C_a^{-1}
///            \end{array} \right]
/// \f]
///
/// The method weighted_model_measure_diff_aug() is another
/// function of the parameters, and the method
/// weighted_jacobian_aug() is the Jacobian of 
/// weighted_model_measure_diff_aug().
///
/// \return The matrix returned by 
///         uncert_weighted_jacobian()
///         augmented at the bottom by the matrix returned by
///         a_priori_cov_chol_inv()
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> weighted_jacobian_aug();


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// The following are tertiary convenient public methods.
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
/// \brief Returns a-posteriori covariance matrix
///
/// After the parameters of the model are retrieved such that the
/// model is fitted to the measurement data, then a call to this 
/// method returns the a-posteriori covariance matrix.
///
/// Assume the following:
///   - K (computed model Jacobian)
///   - Se (measurement error covariance matrix)
///   - Sa (a priori covariance matrix)
///   - J (the matrix returned by weighted_jacobian_aug())
///
/// Then this method returns
/// \f[ 
///     (J^TJ)^{-1} = \left( K^T S_e^{-1} K + S_a^{-1} \right)^{-1}.
/// \f]
///
/// In the context of the Non-Linear Least Squares (NLLS) problem,
/// where J is the Jacobian of the NLLS problem, 
/// \f[ (J^TJ)^{-1} \f]
/// is known as the best fit covariance.
///
/// I am not sure where the best location to implement this 
/// method is.  Maybe it is better to rename it best_fit_covariance()
/// and add to the class NLLSProblem as a method.  Or, maybe there
/// are two different perspectives of the same thing:
///   - a-posteriori covariance
///   - best fit covariance
///
/// If two perspectives, then perhaps it is best to keep 
/// a_posteriori_covariance() method here, and also implement
/// best_fit_covariance() as a member of NLLSProblem class.  If
/// we decide to implement both methods, there will not be a lot of 
/// duplicate code because the body of a_posteriori_covariance()
/// is just a function call.
///
/// \return A-posteriori covariance matrix
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> a_posteriori_covariance();


//-----------------------------------------------------------------------
/// \brief Returns the averaging kernel matrix
///
/// This is another method that I am not sure about the best
/// location for it.  I chose it to be a method of this
/// class because the class has all the data necessary to 
/// compute the averaging kernel.
///
/// Assume the following:
///   - K (computed model Jacobian)
///   - Se (measurement error covariance matrix)
///   - Sa (a priori covariance matrix)
///
/// Then the averaging kernel is
/// \f[
///     \left( K^T S_e^{-1} K + S_a^{-1} \right)^{-1} K^T S_e^{-1} K
/// \f]
///
/// \return The averaging kernel matrix
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> averaging_kernel();


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// More convenient methods for error analysis
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
/// \brief Returns the Cholesky decomposition of the a-priori 
///        covariance matrix
///
/// Assume
///   - Sa (a priori covariance matrix)
///
/// Then the Cholesky decomposition of the a priori covariance matrix is 
/// \f[
///     S_a = C_a C_a^T
/// \f]
/// and this method returns
/// \f[
///     C_a
/// \f]
///
/// \return The Cholesky decomposition of the a-priori 
///         covariance matrix
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> a_priori_cov_chol() const
  { return Sa_chol.copy(); }


//-----------------------------------------------------------------------
/// \brief Returns the square root of the diagonal of the a-priori
///        covariance matrix
///
/// A-priori covariance matrix is returned by a_priori_cov() method,
/// and this method returns the square root of the diagonal of the 
/// matrix.
///
/// \return The square root of the diagonal of the a-priori
///         covariance matrix
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> param_a_priori_uncertainty() const;


//-----------------------------------------------------------------------
/// \brief Returns the square root of the diagonal of the a-posteriori
///        covariance matrix.
///
/// A-posteriori covariance matrix is returned by 
/// a_posteriori_covariance() method, and this method returns the
/// square root of the diagonal of the matrix.
///
/// I am not sure where the best location for implementing this 
/// method is; however, it should be where a_posteriori_covariance()
/// method is.
///
/// \return The square root of the diagonal of the a-posteriori
///         covariance matrix
//-----------------------------------------------------------------------
  virtual blitz::Array<double, 1> param_a_posteriori_uncertainty();


protected:

  blitz::Array<double, 1> Xa;
  blitz::Array<double, 2> Sa;

  // For convenience
  blitz::Array<double, 2> Sa_chol;
  blitz::Array<double, 2> Sa_chol_inv;


};
}
#endif
