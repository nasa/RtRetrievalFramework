#ifndef CHAPMAN_BOA_H
#define CHAPMAN_BOA_H

#include <blitz/array.h>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "auto_derivative.h"
#include "refractive_index.h"

namespace FullPhysics {

void straightline_geometry( // Inputs
			   bool do_plane_parallel,
			   double rearth,
			   const blitz::Array<AutoDerivative<double>, 1>& heights,
			   const blitz::Array<double, 1>& sza_values,
			   // Outputs
			   blitz::Array<AutoDerivative<double>, 2>& chapmanfactors,
			   blitz::Array<AutoDerivative<double>, 1>& toa_nadir_szangles,
			   blitz::Array<AutoDerivative<double>, 1>& toa_entry_szangles);

void refractive_geometry( // Input
			 double rearth,
			 const blitz::Array<AutoDerivative<double>, 1>& heights,
			 const std::vector<boost::shared_ptr<AtmRefractiveIndex> >& refr_index,
			 const blitz::Array<double, 1>& sza_values,
			 // Output
			 blitz::Array<AutoDerivative<double>, 2>& chapmanfactors,
			 blitz::Array<AutoDerivative<double>, 1>& toa_nadir_szangles,
			 blitz::Array<AutoDerivative<double>, 1>& toa_entry_szangles);

void refractive_bending( // Input
			int start_lay, //  Starting at Level n and working out to TOA
			double theta,  //  Function Guess (level angle)
			double dtr,
			const blitz::Array<AutoDerivative<double>, 1>& rtimesn, 
			const blitz::Array<int, 1>& nquad,
			//  layer fine gridding
			const blitz::Array<AutoDerivative<double>, 2>& fineweight, 
			const blitz::Array<AutoDerivative<double>, 2>& fineradii,
			const blitz::Array<AutoDerivative<double>, 2>& finersqnsq, 
			//  Output local distances, cumulative angles and function value
			blitz::Array<AutoDerivative<double>, 1>& distances, 
			AutoDerivative<double>& alpha, 
			AutoDerivative<double>& anglefunc );

/// Output
/// ======
/// partial and whole layer fine gridding:
/// nquad
///
/// Fine layering stuff:
/// fineweight, fineradii, finersqnsq

void full_gridding( // Inputs
		   double rearth,
		   //  Gridding Control
		   int nlower,
		   int nupper,
		   double changeheight,
		   const blitz::Array<AutoDerivative<double>, 1>& heights,
		   const std::vector<boost::shared_ptr<AtmRefractiveIndex> >& refr_index,
		   // Outputs
		   //  Fine layering stuff
		   blitz::Array<AutoDerivative<double>, 2>& fineweight,
		   blitz::Array<AutoDerivative<double>, 2>& fineradii,
		   blitz::Array<AutoDerivative<double>, 3>& finersqnsq,
		   //  partial and whole layer fine gridding
		   blitz::Array<int, 1>& nquad);

/****************************************************************//**
  This class computes Bottom of the Atmosphere radiance
*******************************************************************/

class ChapmanBOA { 

public:

  ChapmanBOA(double rearth,
	     bool do_plane_parallel,
	     const blitz::Array<double, 1>& sza_values,
	     const blitz::Array<AutoDerivative<double>, 1>& heights);

  ChapmanBOA(double rearth,
	     const blitz::Array<double, 1>& sza_values,
	     const blitz::Array<AutoDerivative<double>, 1>& heights,
	     const std::vector<boost::shared_ptr<AtmRefractiveIndex> >& refr_index);

  ChapmanBOA(double rearth,
	     const double sza,
	     const blitz::Array<AutoDerivative<double>, 1>& heights,
	     const boost::shared_ptr<AtmRefractiveIndex>& refr_index);

  // Access angles set up during class construction
  const blitz::Array<AutoDerivative<double>, 2>& chapman_factors() const { return chapmanfactors_; }
  const blitz::Array<AutoDerivative<double>, 1>& toa_nadir_szangles() const { return toa_nadir_szangles_; }
  const blitz::Array<AutoDerivative<double>, 1>& toa_entry_szangles() const  { return toa_entry_szangles_; }

  /// Compute transmittance using molecular extinction and pre-computed angles
  const blitz::Array<AutoDerivative<double>, 1> transmittance(const blitz::Array<AutoDerivative<double>, 1>& extinction) const;
  const AutoDerivative<double> transmittance(const blitz::Array<AutoDerivative<double>, 1>& extinction, int beam_index) const;
    
private:

  blitz::Array<AutoDerivative<double>, 2> chapmanfactors_;
  blitz::Array<AutoDerivative<double>, 1> toa_nadir_szangles_;
  blitz::Array<AutoDerivative<double>, 1> toa_entry_szangles_;

};
}
#endif
