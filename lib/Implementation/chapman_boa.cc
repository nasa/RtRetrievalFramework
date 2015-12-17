#include "chapman_boa.h"
#include "rf_gauleg.h"
#include "refractive_index.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

/// Straight line geometry constructor
ChapmanBOA::ChapmanBOA(double rearth,
		       bool do_plane_parallel,
		       const blitz::Array<double, 1>& sza_values,
		       const blitz::Array<AutoDerivative<double>, 1>& heights)
{
  straightline_geometry(do_plane_parallel, rearth, heights, sza_values, chapmanfactors_, toa_nadir_szangles_, toa_entry_szangles_);
}

/// Refractive geometry constructor
ChapmanBOA::ChapmanBOA(double rearth,
		       const blitz::Array<double, 1>& sza_values,
		       const blitz::Array<AutoDerivative<double>, 1>& heights,
		       const std::vector<boost::shared_ptr<AtmRefractiveIndex> >& refr_index)
{
  refractive_geometry(rearth, heights, refr_index, sza_values, chapmanfactors_, toa_nadir_szangles_, toa_entry_szangles_);
}

/// Refractive geometry but just one set of angles
ChapmanBOA::ChapmanBOA(double rearth,
		       const double sza,
		       const blitz::Array<AutoDerivative<double>, 1>& heights,
		       const boost::shared_ptr<AtmRefractiveIndex>& refr_index) 
{
  Array<double, 1> sza_values(1);
  sza_values(0) = sza;

  std::vector<boost::shared_ptr<AtmRefractiveIndex> > refr_index_arr;
  refr_index_arr.push_back(refr_index);

  refractive_geometry(rearth, heights, refr_index_arr, sza_values, chapmanfactors_, toa_nadir_szangles_, toa_entry_szangles_);
}

/// Compute transmittance from extinction for all configured angles
const blitz::Array<AutoDerivative<double>, 1> ChapmanBOA::transmittance(const blitz::Array<AutoDerivative<double>, 1>& extinction) const
{
  Array<AutoDerivative<double>, 1> trans(chapmanfactors_.extent(secondDim));
  for(int ib = 0; ib < chapmanfactors_.extent(secondDim); ib++) {
    trans(ib) = transmittance(extinction, ib);
  }
  return trans;
}

/// Compute transmittance from extinction for a specific angle index
const AutoDerivative<double> ChapmanBOA::transmittance(const blitz::Array<AutoDerivative<double>, 1>& extinction, int beam_index) const
{
  AutoDerivative<double> trans = 1.0e0;
  for(int n = 0; n < chapmanfactors_.extent(firstDim); n++) {
    AutoDerivative<double> tr  = 0.0e0;
    AutoDerivative<double> arg = extinction(n) * chapmanfactors_(n,beam_index);
    if ( arg < 32.0e0 ) tr = exp ( - arg );
    trans = trans * tr;
  }
  return trans;
}

/****************************************************************//**
 Free functions translated from Fortran
*******************************************************************/

void FullPhysics::straightline_geometry( // Inputs
					bool do_plane_parallel,
					double rearth,
					const blitz::Array<AutoDerivative<double>, 1>& heights,
					const blitz::Array<double, 1>& sza_values,
					// Outputs
					blitz::Array<AutoDerivative<double>, 2>& chapmanfactors,
					blitz::Array<AutoDerivative<double>, 1>& toa_nadir_szangles,
					blitz::Array<AutoDerivative<double>, 1>& toa_entry_szangles ) 
{

  // Sizes of input data dimensions
  int nlayers = heights.extent(firstDim)-1;
  int nbeams = sza_values.extent(firstDim);

  // Resize output arguments
  chapmanfactors.resize(nlayers,nbeams);
  toa_nadir_szangles.resize(nbeams);
  toa_entry_szangles.resize(nbeams);
 
  Array<AutoDerivative<double>, 1> radii(nlayers+1);
      
  // initialize
  toa_nadir_szangles = 0.0e0;
  toa_entry_szangles = 0.0e0;
  chapmanfactors = 0.0e0;

  // set up some local values
  double dtr = atan(1.0e0)/45.0e0;

  //  set up local atmosphere quantities
  for(int n = 0; n < nlayers+1; n++)
    radii(n) = rearth + heights(n);

  //  start Loop over solar beams
  //  ===========================
  for(int ibeam = 0; ibeam < nbeams; ibeam++) {

    //  get TOA solar zenith angle
    double xboa = sza_values(ibeam);
    toa_nadir_szangles(ibeam) = xboa;
    double xboa_r = xboa * dtr;
    double mu_boa = cos(xboa_r);

    if ( do_plane_parallel ) {
      //  Plane parallel

      toa_entry_szangles (ibeam) = xboa;
      for(int k = 0; k < nlayers; k++)
	chapmanfactors(k,ibeam) = 1.0e0 / mu_boa;

    } else {
      //  Curved atmopshere
      //    sine-rule; PHI = earth-centered angle

      double gm_boa = sqrt(1.0e0 - mu_boa*mu_boa);
      AutoDerivative<double> sinth1 = gm_boa * radii(nlayers) / radii(0);
      AutoDerivative<double> sth1   = asin(sinth1);
      toa_entry_szangles(ibeam) = sth1 / dtr;
      AutoDerivative<double> re_upper = radii(0);

      for(int k = 0; k < nlayers; k++) {
	AutoDerivative<double> delz = heights(k) - heights(k+1);
	AutoDerivative<double> re_lower = re_upper - delz;
	AutoDerivative<double> sinth2 = re_upper * sinth1 / re_lower;
	AutoDerivative<double> sth2   = asin(sinth2);
	AutoDerivative<double> phi    = sth2 - sth1;
	AutoDerivative<double> sinphi = sin(phi);
	AutoDerivative<double> dist = re_upper * sinphi / sinth2;
	chapmanfactors(k,ibeam) = dist / delz;

	re_upper = re_lower;
	sinth1 = sinth2;
	sth1   = sth2;
      }
    }

  }    //  End beam loop
   
  // end subroutine straightline_geometry
}

void FullPhysics::refractive_geometry( // Input
				      double rearth,
				      const blitz::Array<AutoDerivative<double>, 1>& heights,
				      const std::vector<boost::shared_ptr<AtmRefractiveIndex> >& refr_index,
				      const blitz::Array<double, 1>& sza_values,
				      // Output
				      blitz::Array<AutoDerivative<double>, 2>& chapmanfactors,
				      blitz::Array<AutoDerivative<double>, 1>& toa_nadir_szangles,
				      blitz::Array<AutoDerivative<double>, 1>& toa_entry_szangles)
{
  // Sizes of input data dimensions
  int nlayers = heights.extent(firstDim)-1;
  int nbeams = sza_values.extent(firstDim);

  if((int) refr_index.size() != nbeams) {
    Exception err;
    err << "Number of refractive index objects: " << refr_index.size()
	<< " does not match number of solar zenith angles: " << nbeams;
    throw err;
  }
     
  // Resize output arguments

  //  chapman factors to BOA, and TOA Sza angles
  chapmanfactors.resize(nlayers,nbeams);
  toa_nadir_szangles.resize(nbeams);
  toa_entry_szangles.resize(nbeams);

  //  Local
  //  =====

  //  Set initial gridding Control variables
  int nupper = 10;
  int nlower = 20;
  double changeheight = 20.0e0;
  
  // Will never allocate more than this
  int max_nquads = 100;//max(nupper, nlower);

  //  ray constants
  Array<AutoDerivative<double>, 2> rtimesn(nlayers+1,nbeams);
  Array<AutoDerivative<double>, 1> radii(nlayers+1);
 
  //  partial and whole layer fine gridding
  Array<int, 1> nquad(nlayers);

  //  Fine layering stuff (whole atmosphere)
  Array<AutoDerivative<double>, 2> fineweight (nlayers,max_nquads);
  Array<AutoDerivative<double>, 2> fineradii  (nlayers,max_nquads);
  Array<AutoDerivative<double>, 3> finersqnsq (nlayers,max_nquads,nbeams);

  //  initialize
  toa_nadir_szangles = 0.0e0;
  toa_entry_szangles = 0.0e0;
  chapmanfactors     = 0.0e0;
  Range ra = Range::all();

  //  set up some local values
  double dtr = atan(1.0e0)/45.0e0;

  //  set up local atmosphere quantities
  for (int i = 0; i <= nlayers; i++) {
    radii(i)   = rearth + heights(i);
    for(int ibeam = 0; ibeam < nbeams; ibeam++)
      rtimesn(i, ibeam) = radii(i) * refr_index[ibeam]->at_level(i);
  }

  //  full gridding of atmosphere
  full_gridding                                         
    ( rearth, nlower, nupper, changeheight, heights, refr_index, 
      fineweight, fineradii, finersqnsq, nquad );

  //  start Loop over solar beams
  //  ===========================
  for(int ibeam = 0; ibeam < nbeams; ibeam++) {

    //  local solar zenith angles
    double xboa = sza_values(ibeam);

    //  refraction from BOA, with distances
    Array<AutoDerivative<double>, 1> distances(nlayers);
    AutoDerivative<double> xtoa, alpha; // returned values
    refractive_bending
      ( nlayers, xboa, dtr, rtimesn(ra, ibeam), nquad, 
	fineweight, fineradii, finersqnsq(ra, ra, ibeam),
	distances, alpha, xtoa );

      //  get TOA solar zenith angle, BOA-level Chapmans
    toa_nadir_szangles (ibeam) = xtoa;
    toa_entry_szangles (ibeam) = alpha;
    for(int k = nlayers-1; k >= 0; k--) {
      AutoDerivative<double> delz = heights(k) - heights(k+1);
      chapmanfactors(k,ibeam) = distances(k) / delz;
    }

  } //  End beam loop

} //  End of refractive_geometry

void FullPhysics::refractive_bending( // Input
				     int start_lay, //  Starting at Level n and working out to TOA
				     double theta,  //  Function Guess (level angle)
				     double dtr,
				     const blitz::Array<AutoDerivative<double>, 1>& rtimesn, 
				     const blitz::Array<int, 1>& nquad,
				     //  Layer fine gridding
				     const blitz::Array<AutoDerivative<double>, 2>& fineweight, 
				     const blitz::Array<AutoDerivative<double>, 2>& fineradii,
				     const blitz::Array<AutoDerivative<double>, 2>& finersqnsq, 
				     //  Output local distances, cumulative angles and function value
				     blitz::Array<AutoDerivative<double>, 1>& distances, 
				     AutoDerivative<double>& alpha, 
				     AutoDerivative<double>& anglefunc )
{

  // Sizes of input data dimensions
  int nlayers = rtimesn.extent(firstDim)-1;

  // Resize output arrays
  distances.resize(nlayers);

  //  Initialize
  anglefunc = 0.0e0;
  distances = 0.0e0;

  //  constants
  double theta_r  = theta*dtr;
  double sintheta = sin(theta_r);
  AutoDerivative<double> raycons = rtimesn(start_lay)*sintheta;
  AutoDerivative<double> csq     = raycons * raycons;
  alpha = asin(raycons/rtimesn(0))/dtr;
  
  //  distance calculation
  for(int i = start_lay-1; i >= 0; i--) {
    AutoDerivative<double> sum = 0.0e0;
    for(int k = 0; k < nquad(i); k++) {
      AutoDerivative<double> Fik     = 1.0e0 / ( finersqnsq(i,k) - csq);
      AutoDerivative<double> Xik_sq  = 1.0e0 + csq * Fik;
      AutoDerivative<double> Xik     = sqrt(Xik_sq);
      sum = sum + fineweight(i,k) * Xik;
      distances(i) = sum;
    }
  }

  //  angle calculations
  AutoDerivative<double> phi_rad = 0.0e0;
  for(int i = start_lay-1; i >= 0; i--) {
    AutoDerivative<double> sum = 0.0e0;
    for(int k = 0; k < nquad(i); k++) {
      AutoDerivative<double> help = sqrt((finersqnsq(i,k)/csq) - 1.0);
      sum = sum + fineweight(i,k)/help/fineradii(i,k);
    }
    phi_rad = phi_rad + sum;
  }
  AutoDerivative<double> phicum = phi_rad / dtr;
  
  //  SZA at beam entry + total E-C angle should equal SZA at TOA
  //   This function should  be zero
  //      anglefunc  = theta + phicum - alpha
  anglefunc  = phicum + alpha;

} // End of refractive_bending

void FullPhysics::full_gridding( // Inputs
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
				blitz::Array<int, 1>& nquad)
{
  
  //  Input arguments
  //  ---------------

  //  Dimensioning
  int nlayers = heights.extent(firstDim)-1;
  int max_nquads = fineweight.extent(secondDim);

  //  initialize
  nquad = 0; 
  fineweight = 0.0e0;
  fineradii = 0.0e0; 
  finersqnsq = 0.0e0;
  
  //  Set up gridding
  for(int i = 0; i < nlayers; i++) {
    nquad(i) = nlower;
    if (heights(i) > changeheight)
      nquad(i) = nupper;
  }

  Array<AutoDerivative<double>, 1> xgrid(max_nquads);
  Array<AutoDerivative<double>, 1> wgrid(max_nquads);

  //  Gauss-Legendre  Integration
  for(int i = 0; i < nlayers; i++) {
    rf_gauleg(heights(i+1), heights(i), xgrid, wgrid, nquad(i));

    for(int k = 0; k < nquad(i); k++) {
      fineradii(i,k)  = rearth + xgrid(k);
      fineweight(i,k) = wgrid(k);
      AutoDerivative<double> xd = (xgrid(k) - heights(i+1)) / (heights(i) - heights(i+1));
      for(int beam_idx = 0; beam_idx < finersqnsq.extent(thirdDim); beam_idx++) {
	AutoDerivative<double> fn = refr_index[beam_idx]->at_layer(i, xd);
	finersqnsq(i,k,beam_idx) = pow(fineradii(i,k), 2) * pow(fn, 2);
      }
    }
  }
} // End of full_gridding
