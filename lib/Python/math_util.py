from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
import math
import bisect
import copy
import six

import numpy

def smooth(x,window_len=10,window='hanning'):
    """smooth the data using a window with requested size.

    Downloaded from: http://www.scipy.org/Cookbook/SignalSmooth
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]

    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval(window+'(window_len)')

    y=convolve(old_div(w,w.sum()),s,mode='same')
    return y[window_len-1:-window_len+1]

def linear_interpol(data_in, grid_in, grid_out, extrapolate=False, extrap_err=True):
    'Perform a linear interpolation from one grid to another'

    if len(grid_in) == 0:
        raise ValueError('grid_in can not be empty')

    if len(grid_out) == 0:
        raise ValueError('grid_out can not be empty')

    new_data = numpy.zeros(len(grid_out), dtype=float)

    # We don't want to extrapolate
    if (grid_in[0] > grid_out[0]) and not extrapolate and extrap_err:
        raise ValueError("Can not interpolate since lowest input grid point: %f is larger than lowest output grid point: %f" % (grid_in[0], grid_out[0]))

    # Look for possible precision problem when bottom levels are near each other
    if (grid_in[len(grid_in) - 1] < grid_out[len(grid_out) - 1]):
        if (abs(grid_in[len(grid_in) - 1] - grid_out[len(grid_out) - 1]) < 1e-16):
            # Deal with percision problem by setting new grid max to old grid max
            # since it seems that they are essentially equal
            grid_out[len(grid_out) - 1] = grid_in[len(grid_in) - 1]
        elif not extrapolate and extrap_err:
            raise ValueError("highest input grid point %f is less than highest output grid grid point %f" % \
                  (grid_in[len(grid_in) - 1], grid_out[len(grid_out) - 1]))

    for out_g_row in range(len(grid_out)):
        y = grid_out[out_g_row]
        row = bisect.bisect_left(grid_in, y)

        row = min([row, len(data_in)-1])
        row = max([row, 0])

        if (row == 0):
            if extrapolate:
                frac = old_div((y - grid_in[row]),(grid_in[row+1] - grid_in[row]))
                new_data[out_g_row] = data_in[row] + frac * \
                                      (data_in[row+1] - data_in[row])
            else:
                new_data[out_g_row] = float('nan')#data_in[row]
        elif (data_in[row] == data_in[row - 1]):
            new_data[out_g_row] = data_in[row]
        else:
            frac = old_div((y - grid_in[row-1]),(grid_in[row] - grid_in[row-1]))
            new_data[out_g_row] = data_in[row - 1] + frac * \
                                  (data_in[row] - data_in[row - 1])

    return new_data

def linear_interpol_idl(v, x, u = None):
    """
    Linearly interpolate vectors with a regular or irregular grid.

    CALLING SEQUENCE:
        result = linear_interpol(V, N)          ;For regular grids.

	result = linear_interpol(V, X, U)       ;For irregular grids.

    INPUTS:
	V:	The input vector can be any type except string.

	For regular grids:
	N:	The number of points in the result when both input and
		output grids are regular.

	Irregular grids:
	X:	The absicissae values for V.  This vector must have same # of
		elements as V.  The values MUST be monotonically ascending
		or descending.

	U:	The absicissae values for the result.  The result will have
		the same number of elements as U.  U does not need to be
		monotonic.  If U is outside the range of X, then the
		closest two endpoints of (X,V) are linearly extrapolated.
    OUTPUTS:
	INTERPOL returns a floating-point vector of N points determined
	by interpolating the input vector by the specified method.

    PROCEDURE:
	For linear interpolation,
	Result(i) = V(x) + (x - FIX(x)) * (V(x+1) - V(x))

	where 	x = i*(m-1)/(N-1) for regular grids.
		m = # of elements in V, i=0 to N-1.

	For irregular grids, x = U(i).
		m = number of points of input vector.

    SOURCE:
         IDL interpol.pro routine
    """

    m = len(v)               # of input pnts

    regular = u is None
    
    if regular:              # Simple regular case?
        r = arange(x, dtype=float)*(old_div((m-1.0),((x-1.0) > 1.0))) #Grid points in V
        rl = int(r)	      # Cvt to integer
        dif = v[1:]-v        # Other types are already signed
            
        return V[rl] + (r-rl)*dif[rl] # interpolate

    if len(x) != m:
        raise ValueError('V and X arrays must have same # of elements')

    p = numpy.zeros(len(u), dtype=float)
    for out_row in range(len(u)):
        s = bisect.bisect_left(x, u[out_row])

        s = max([0, s])
        s = min([s, m-2])

        # Fix the value bisect returns to agree with how IDL's VALUE_LOCATE
        # finds values
        if x[s] > u[out_row]:
            s -= 1

        p[out_row] = (u[out_row]-x[s])*(v[s+1]-v[s])/(x[s+1] - x[s]) + v[s]

    return p

def resample_profile(in_pres, in_data, resample_to, extrapolate=False, log_data=False):

    pmin = in_pres[0]
    pmax = in_pres[len(in_pres)-1]

    in_use_data = numpy.zeros(len(in_data), dtype=float)
    for row in range(len(in_data)):
        loc_value = in_data[row]

        if log_data:
            if abs(loc_value) <= 1.0e-30:
                loc_value = 1.0e-30
            in_use_data[row] = math.log10( loc_value )
        else:
            in_use_data[row] = loc_value
   
    # Create log of source pressure grid
    src_grid_log = numpy.zeros(len(in_pres), dtype=float)
    for i in range (len(in_pres)):
        src_grid_log[i] = math.log10(in_pres[i])
        
    # Create log of destination pressure grid
    if isinstance(resample_to, six.integer_types):
        dst_grid_log = numpy.zeros(resample_to, dtype=float)
        dst_data = numpy.zeros(resample_to, dtype=float)
        for i in range (resample_to):
            dst_grid_log[i] = math.log10(pmin + i * (pmax - pmin) / (resample_to - 1))
    else:
        dst_grid_log = numpy.zeros(len(resample_to), dtype=float)
        dst_data = numpy.zeros(len(resample_to), dtype=float)
        for i in range (len(resample_to)):
            dst_grid_log[i] = math.log10(resample_to[i])

    if extrapolate == True:
        dst_data = linear_interpol(in_use_data, src_grid_log, dst_grid_log, extrapolate=True)
    else:
        top_level = 0
        while dst_grid_log[top_level] < src_grid_log[0] and top_level < len(dst_grid_log)-1:
            top_level += 1
        
        btm_level = 0
        while dst_grid_log[btm_level] < src_grid_log[len(src_grid_log)-1] and btm_level < len(dst_grid_log)-1:
            btm_level += 1
        btm_level -= 1
        
        dst_data[top_level:btm_level+1] = linear_interpol(in_use_data, src_grid_log, dst_grid_log[top_level:btm_level+1], extrapolate=False)

        for level in range(0, top_level+1):
            dst_data[level] = dst_data[top_level+1]

        for level in range(btm_level+1, len(dst_grid_log)):
            dst_data[level] = dst_data[btm_level]
        
    if log_data:
        for row in range(len(dst_data)):
            dst_data[row] = math.pow(10, dst_data[row])           
  
    return dst_data

def resample_aerosol_old(pressures_in, aerosol_data, pressures_out, debug=False):
    # Lowest extinction level
    ext_0 = 1.0E-20

    if type(aerosol_data) == numpy.ndarray and len(aerosol_data.shape) == 2:
        aer_matrix = aerosol_data
    else:
        aer_matrix = numpy.zeros((len(aerosol_data), 1), dtype=float)
        aer_matrix[:,0] = aerosol_data[:]

    # Make sure there are no stray negatives in there..
    aer_matrix = abs(aer_matrix)

    if debug:
        print('input aod = %f, %d levels' % (total_aod(pressures_in, aer_matrix), aer_matrix.shape[0]))
        print('')

    num_in_levels = len(pressures_in)
    num_out_levels = len(pressures_out)
    num_aerosols = aer_matrix.shape[1]

    # Smoothing width -- seems to be tied to gridding of regular grid
    width = num_in_levels * 10

    # create regular press-levels
    num_regular_levels = old_div(int(pressures_in[num_in_levels-1]-pressures_in[0]),10)+2
    pressures_regular = arange(1, num_regular_levels*10, 10, dtype=float)

    # Create output data matrix
    ext_out = numpy.zeros((num_aerosols, num_out_levels), dtype=float)

    for aero_idx in range(num_aerosols):

        # interpolate to regular 10Pa levels
        desired_aod = total_aod(pressures_in, aer_matrix[:,aero_idx])
        aero_interp = linear_interpol_idl(aer_matrix[:,aero_idx], pressures_in, pressures_regular)

        # smooth aerosol profile
        aero_smooth = smooth(aero_interp, width, window='flat')
        #aero_smooth = aero_interp

        if debug:
            print('  type input aod [%d] = %f' % (aero_idx, desired_aod))
            print(' type smooth aod [%d] = %f' % (aero_idx, total_aod(pressures_regular, aero_smooth)))

        btm_level = 0
        while pressures_out[btm_level] < pressures_in[num_in_levels-1] and btm_level < num_out_levels-1:
            btm_level += 1

        # set lowest value to surface pressure
        p1 = numpy.zeros(num_out_levels, dtype=float)
        if type(pressures_out) == numpy.ndarray:
            p1 = pressures_out.copy()
        else:
            p1 = copy(pressures_out)
        p1[btm_level] = pressures_in[num_in_levels-1]

        pressures_fine = numpy.zeros(101, dtype=float)
        ext = numpy.zeros(num_out_levels, dtype=float)
        ext[:] = ext_0

        aer_tot_aod = 0.0
        for out_lev_idx in range(1,btm_level+1):
            del_p = p1[out_lev_idx] - p1[out_lev_idx-1]
            
            for k in range(100+1):
                pressures_fine[k] = p1[out_lev_idx - 1] + (old_div(del_p, 100.0E0)) * float(k)
                
            aero_interp = linear_interpol_idl(aero_smooth, pressures_regular, pressures_fine)
            tot_od = 0.0E0
            
            # integrate optical depth for on layer
            for k in range(100):
                del_od = (old_div(del_p, 100.0E0)) * (aero_interp[k] + aero_interp[k+1]) / 2.0
                tot_od += del_od

            aer_tot_aod += tot_od

            #print ' tot_aod[%d, %d] = %f' % (aero_idx, out_lev_idx, tot_od), ':', del_p
   
            # compute extinction of layer level
            ext[out_lev_idx]= 2.0 * tot_od / del_p - ext[out_lev_idx-1]
            if ext[out_lev_idx] < ext_0:
                ext[out_lev_idx] = ext_0

        ext_out[aero_idx,:] = linear_interpol_idl(ext, p1, pressures_out)

        new_aod = total_aod(pressures_out, ext_out[aero_idx,:])       
        ext_scale = old_div(desired_aod,new_aod)

        if debug:
            print(' type interp aod [%d] = %f * %f scaling' % (aero_idx, new_aod, ext_scale))
        
        ext_out[aero_idx,:] = ext_out[aero_idx,:] * ext_scale

        for lev_idx in range(num_out_levels):
            if ( ext_out[aero_idx, lev_idx] <= ext_0):
                ext_out[aero_idx, lev_idx] = ext_0

        if debug:
            print(' type scaled aod [%d] = %f' % (aero_idx, total_aod(pressures_out, ext_out[aero_idx,:])))
            print('')

    if debug:
        print(' destination aod = %f, %d levels' % (total_aod(pressures_out, ext_out.transpose()), ext_out.shape[1]))

    return ext_out

def resample_aerosol(plev1, alev1_in, plev2, log_data=False, tau1=None, tau2=None, debug=False):
    # NOTE: MUST MODIFY TO ACT DIRECTORY ON AEROSOL AND PRESSURE FILES.
    # WE MAY NEED TO ELIMINATE SOME LAYERS ON OUTPUT IF THEY DO NOT EXIST
    # Thus, PLEV2 must be modified to have the correct surface pressure.
    #
    # INPUTS
    # Plev1: Pressure set 1, on levels [Pa]  (n1)
    # Plev2: Pressure set 2, on levels [Pa]  (n2)
    # Alev1: Aerosol concentrations, set 1, on levels [per Pa]. ; (n1, ntypes)
    #
    # RETURN VALUE
    # Alev2: Aerosol concentrations, set 2, on levels [per Pa]  ; (n2, ntypes)
    #
    # OPTIONAL OUTPUTS
    # Tau1 :
    #
    # This routine will do an interpolation of aerosol concentration
    # from one pressure set to another.
    #
    # It is assumed that the two sets have the same surface pressure

    # initializations
    n1 = len(plev1)
    n2 = len(plev2)
    ntypes = alev1_in.shape[1]
    tau1 = numpy.zeros((n1-1,ntypes), dtype=float)
    tau2 = numpy.zeros((n2-1,ntypes), dtype=float)
    alay2 = numpy.zeros((n2-1,ntypes), dtype=float)
    alev2 = numpy.zeros((n2,ntypes), dtype=float)

    # convert to linear space if necessary
    if log_data:
        alev1 = exp(alev1_in)
    else:
        alev1 = alev1_in

    # STEP 1: Calculate Tau1 (on layers)
    alay1 = 0.5*(alev1[0:n1-1,:] + alev1[1:,:])
    dp1 = plev1[1:] - plev1[0:n1-1]
    for aer_idx in range(0,ntypes):
        tau1[:,aer_idx] = dp1 * alay1[:, aer_idx]
        
    # STEP2 : Distribute Tau1's on new pressure grid, dividing up overlapping layers
    f = linear_interpol_idl(numpy.arange(n1, dtype=float), plev1, plev2)
    j = 0 # index of which hi-res layer which are distributing
    for i in range(0, n2-1): # cycle through output LAYERS
        k1 = int(f[i]) # get rounded-down version of f[i]
        k2 = int(f[i+1])
        nk = k2-k1+1
        for j in range(0,nk): # ex: j=2,3,4,5
            k = j + k1
            first = max([k, f[i]])
            next  = min([k+1, f[i+1]])
            frac  = next-first
            # make sure this layer exists
            if (k < 0 or k > (n1-2)):
                continue
            # add contribution of tau from this layer
            tau2[i,:] += tau1[k,:] * frac

    # STEP3 : Determine new aerosol concentration on Layers
    dp2 = plev2[1:] - plev2[0:n2-1]
    for a in range(0, ntypes):
        alay2[:,a] = old_div(tau2[:,a], dp2)

    # STEP4 : Interpolate layers -> levels
    # make boundary values the same, average in between
    alev2[0,:]      = alay2[0,:]
    alev2[n2-1,:]   = alay2[n2-2,:]
    alev2[1:n2-1,:] = 0.5 * (alay2[0:n2-2,:] + alay2[1:n2-1,:])

    # convert to log space if necessary
    if log_data:
        alev2 = alog(alev2 > 1e-20)

    # Put back into indexing used by python
    return numpy.transpose(alev2)

def total_aod(pressure, aerosol, return_layers=False, in_log=False):

  if type(pressure) == numpy.ndarray:
      num_levels = pressure.shape[0]
  else:
      num_levels = len(pressure)

  if type(aerosol) == numpy.ndarray and len(aerosol.shape) == 2:
      aer_matrix = aerosol
  else:
      aer_matrix = numpy.zeros((len(aerosol), 1), dtype=float)
      aer_matrix[:,0] = aerosol[:]

  # convert from log if aerosols are in log. 
  # If there are any negatives (should all since no aerosol > 1
  # then we have log retrieval of aerosol
  if in_log:
      aer_matrix = numpy.exp(aer_matrix)

  if aer_matrix.shape[0] != num_levels:
      raise ValueError('Number of pressure levels does not match number of aerosol levels')

  num_aer_types = aer_matrix.shape[1]
      
  # compute total aod
  layer_aod   = numpy.zeros(num_levels-1, dtype=float)
  layer_press = numpy.zeros(num_levels-1, dtype=float)
  for layer_idx in range(num_levels-1):
      delta_press = old_div((pressure[layer_idx + 1] - pressure[layer_idx]), 2.0)
      layer_press[layer_idx] = pressure[layer_idx + 1] - delta_press
      for aer_idx in range(num_aer_types):
          layer_aod[layer_idx] = layer_aod[layer_idx] + \
                                 (delta_press *
                                  (aer_matrix[layer_idx, aer_idx] + aer_matrix[layer_idx + 1, aer_idx]))

  total_aod = sum(layer_aod)

  if return_layers:
      return layer_aod, layer_press
  else:
      return total_aod

def jpl_gravity(gdlat, altit):
    """Computes the effective Earth gravity at a given latitude and altitude.
    This is the sum of the gravitational and centripital accelerations.
    These are based on equation I.2.4-(17) in US Standard Atmosphere 1962
    The Earth is assumed to be an oblate ellipsoid, with a ratio of the
    major to minor axes = sqrt(1+con) where con=.006738
    This eccentricity makes the Earth's gravititational field smaller at the
    poles and larger at the equator than if the Earth were a sphere of the
    same mass. It also makes the local mid-latitude gravity field not point
    toward the center of mass.

    Input Parameters:
        gdlat       r*4  GeoDetric Latitude (degrees)
        altit       r*4  Geometric Altitude (km)

    Output Parameter:
        gravity     r*4  Effective Gravitational Acceleration (m/s2)

    Interestingly, since the centripital effect of the Earth's rotation
    (-ve at equator, 0 at poles) has almost the opposite shape to the
    second order gravitational field (+ve at equator, -ve at poles), their
    sum is almost constant so that the surface gravity can be approximated
    (to .07%) by the simple expression g = 0.99746*GM/radius**2, the latitude
    variation coming entirely from the variation of surface r with latitude.
    """
    d2r=old_div(3.14159265,180.)
    gm=3.9862216e+14
    omega=7.292116E-05
    con=.006738
    shc=1.6235e-03
    eqrad=6378178.

    # Convert from geodetic latitude (GDLAT) to geocentric latitude (GCLAT).
    gclat=numpy.arctan(old_div(numpy.tan(d2r*gdlat),(1.0+con)))  # radians
    # On computers which crash at the poles try the following expression
    # gclat=d2r*gdlat-con*sin(d2r*gdlat)*cos(d2r*gdlat)/(1+con*cos(d2r*gdlat)**2)
    radius=1000.0*altit+old_div(eqrad,numpy.sqrt(1.+con*numpy.sin(gclat)**2))
    ff=(old_div(radius,eqrad))**2
    hh=radius*omega**2
    ge_=old_div(gm,eqrad**2)                       # = gravity at Re
    gravity=(ge_*(1-shc*(3*numpy.sin(gclat)**2-1)/ff)/ff-hh*numpy.cos(gclat)**2) *(1+0.5*(numpy.sin(gclat)*numpy.cos(gclat)*(old_div(hh,ge_)+2*shc/ff**2))**2)

    return gravity

def gravity_profile(P, T, q, elevation, lat):
    """p is ordered space to surface, in Pascals
    elevation is in km
    T is temperature [K] on levels
    q is specific humidity on levels
    lat is in degrees
    
    returns the gravity profile [m/s^2] on LAYERS
    """
    efact = 0.60776532
    Rd = 287.04742
    nz = len(P) - 1
    z = P*0.
    z[nz] = elevation
    grav = numpy.zeros(nz, dtype=float)
    for i in range(nz-1, -1, -1):
        qbar = (q[i] + q[i+1])*0.5
        Pbar = (P[i] + P[i+1])*0.5
        Tbar = (T[i] + T[i+1])*0.5
        zlower = z[i+1]
        dP = P[i+1] - P[i]
        Tv = Tbar * (1.0e0 + qbar * efact)
        logratio = numpy.log(old_div(P[i+1], P[i]))
        grav[i] = jpl_gravity(lat, zlower)
        dz = logratio * Tv * Rd / grav[i] * 1e-3 # km
        grav[i] = jpl_gravity(lat, zlower + 0.5*dz)
        dz = logratio * Tv * Rd / grav[i] * 1e-3 # km
        z[i] = zlower + dz

    return grav

def compute_xco2(co2, p, t=None, h2o=None, ak=None):
    """COMPUTES XCO2 from Profiles of CO2 and P on levels [TOA->surface ordering]
    If T, H2O are present, the water vapor correction will be made.
    (gravity must be computed in this case)
    
    Can take into account optional averaging kernel (the result will be unnormalized!)
    AK is the averaging kernel on levels

    co2 : profile of co2 on levels
    p   : pressure profile [Pa] on levels
    t   : temp profile [K] on levels
    h2o : h2o volume concentration profile on levels"""

    nlev = len(co2)

    # Get pressure weighting function on layers
    if len(p.shape) == 2:
        p = p[:, 0]

    if len(co2.shape) == 2:
        co2 = co2[:, 0]

    # Compute h2o correction if desired
    if t != None and h2o != None:
        if len(t.shape) == 2:
            t = t[:, 0]

        if len(h2o.shape) == 2:
            h2o = h2o[:, 0]

        epsilon = 0.622
        q = epsilon * h2o / (1. + epsilon * h2o)
        qbar = (q[1:] + q[0:-1]) * 0.5
        gbar = gravity_profile(p, t, q, 0., 45.)
        corr = old_div((1-qbar),gbar)
    else:
        corr = numpy.ones(nlev-1)

    dp = p[1:] - p[0:-1]

    wgt = dp * corr
    wgt = old_div(wgt, numpy.sum(wgt))

    # Get the pressure weighting function on levels
    wgtlev = numpy.zeros(nlev, dtype=float)

    wgtlev[0]        = 0.5 * wgt[0]
    wgtlev[1:nlev-1] = 0.5 * (wgt[1:nlev-1] + wgt[0:nlev-2] )
    wgtlev[nlev-1]   = 0.5 * wgt[nlev-2]

    # compute XCO2
    if ak == None:
        ak = numpy.ones(nlev, dtype=float)
    elif len(ak.shape) == 2:
        ak = ak[:,0]

    xco2 = numpy.sum(wgtlev * co2 * ak)

    return xco2




