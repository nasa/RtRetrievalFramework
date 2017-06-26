#!/usr/bin/env python

from __future__ import division
from builtins import zip
from builtins import str
from builtins import range
from builtins import object

import logging
from optparse import OptionParser
from glob import glob
from collections import namedtuple
from datetime import datetime, timedelta

from netCDF4 import Dataset as NetCDFFile
from numpy import *

logger = logging.getLogger()

# Constants
# Equatorial radius (km)
EARTH_RADIUS = 6378.137 # (km)
U_GAS_CONST = 8.31432
AVAGADRO = 6.02217e+23

# Array of altitudes onto which model data are to be interpolated
Z_GRID = arange(0,71)

# Mean molecular weights
MEAN_MOL_WGT = zeros(Z_GRID.shape[0]) + 28.964

# Named tuple types:

# Create named tuple type used for collecting ncep data together
NCEPData = namedtuple("NCEPData", ("air_temp", "spec_hum", "trop_pres", "geo_height", "surface_pres"))
ModelData = namedtuple("ModelData", ["name", "units", "data"])
Model = namedtuple("Model", ["pressure", "temperature", "height", "h2o_vmr", "surface_pressure"])

#*************************************************************
# Interpolates in time, lat/longitude into the NCEP re-analyses to
# produce a model file for a particular site for a series of dates.
#
# Based on the IDL routine mod_maker8.pro and gfit routines
# by Geoff Toon
# This code should match what is produced by gsetup
# from the April 2012 version of GGG:
# https://tccon-wiki.caltech.edu/Software/GGG/Download/GGG_2012_April_Update
#
# Requires that the user has acquired the four necessary .nc files
# (Air Temperature, Geopotential Height, Tropopause Pressure, Specific Humidity)
# from the following four NOAA website:
#
# http://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?&Dataset=NCEP+Reanalysis+Pressure+Level&Dataset=NCEP+Reanalysis+Daily+Averages+&Pressure+Level&Variable=Air+Temperature
# http://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?&Dataset=NCEP+Reanalysis+Pressure+Level&Dataset=NCEP+Reanalysis+Daily+Averages+&Pressure+Level&Variable=Specific+Humidity
# http://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NCEP+Reanalysis+Tropopause+Level&Dataset=NCEP+Reanalysis+Daily+Averages+Tropopause+Level&Variable=Pressure
# http://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?&Dataset=NCEP+Reanalysis+Pressure+Level&Dataset=NCEP+Reanalysis+Daily+Averages+Pressure+Level&Variable=Geopotential+height

# Subroutine svp_wv_over_ice
# Uses the Goff-Gratch equation to calculate the saturation vapor
# pressure of water vapor over ice at a user-specified temperature.
# Input:  temp (K)
# Output: svp (mbar)
def svp_wv_over_ice(temp):
  t0=273.16   # triple point temperature (K)
  tr=t0/temp
  yy=-9.09718*(tr-1)-3.56654*log10(tr)+0.876793*(1-1/tr)
  svp=6.1173*10**yy   # saturation vapor pressure over ice (mbar)
  return svp

# Define US Standard Atmosphere (USSA) for use above 10 mbar
p_ussa=[10.0,  5.0,   2.0,   1.0,   0.1,   0.01,  0.001, 0.0001]
t_ussa=[227.7, 239.2, 257.9, 270.6, 231.6, 198.0, 189.8, 235.0]
z_ussa=[31.1,  36.8,  42.4,  47.8,  64.9,  79.3,  92.0,  106.3]

# SUBROUTINE LINEAR_INTERP0
# Evaluates  fout = fin(xx) for fout of dimension 0 (scalar)
def linear_interp0(fin, xx):
   index_xx = int(xx)
   if index_xx < 0:
       index_xx = 0
   fr_xx=xx-index_xx
   fout=(1-fr_xx)*fin[index_xx]+fr_xx*fin[index_xx+1]
   return fout

# SUBROUTINE LINEAR_INTERP1
# Evaluates  fout = fin(xx) for fout of dimension 1 (vector)
def linear_interp1(fin, xx):
   index_xx = int(xx)
   if index_xx < 0:
       index_xx = 0
   fr_xx=xx-index_xx
   fout=(1-fr_xx)*fin[index_xx,:]+fr_xx*fin[index_xx+1,:]
   return fout

# SUBROUTINE BILINEAR_INTERP1
# Evaluates  fout = fin(xx,yy) for fout of dimension 1 (vector)
def bilinear_interp1(fin, xx, yy):
  index_xx=int(xx)
  fr_xx=xx-index_xx
  index_yy=int(yy)
  fr_yy=yy-index_yy
  fout = \
    (fin[:,index_yy,index_xx] * (1.0 - fr_xx) + \
     fin[:,index_yy,index_xx+1] * fr_xx) * (1.0 - fr_yy) + \
    (fin[:,index_yy+1,index_xx] * (1.0 - fr_xx) + \
     fin[:,index_yy+1,index_xx+1] * fr_xx) * fr_yy
  return fout

# SUBROUTINE BILINEAR_INTERP2
# Evaluates  fout = fin(xx,yy) for fout of dimension 2 (array)
def bilinear_interp2(fin, xx, yy):
  index_xx=int(xx)
  fr_xx=xx-index_xx
  index_yy=int(yy)
  fr_yy=yy-index_yy
  fout= \
    (fin[:,:,index_yy,index_xx] * (1.0 - fr_xx) + \
     fin[:,:,index_yy,index_xx+1] * fr_xx) * (1.0 - fr_yy) + \
    (fin[:,:,index_yy+1,index_xx] * (1.0 - fr_xx) + \
     fin[:,:,index_yy+1,index_xx+1] * fr_xx) * fr_yy
  return fout

def log1pxox(x):
    """Computes log(1+x)/x without risk of underflow or zero-division. For x < 0.1
       LOG1PXOX is approximated by the series 1-x/2+x**2/3-x**3/4+x**4/5-x**5/6"""

    if (abs(x) > 0.1):
        log1pxox = log(1 + x) / x
    else:
        nterm=13  # double precision
        #nterm = 6   # single precision
        log1pxox = 1.e0 / nterm
        #for j = nterm-1,1,-1:
        for j in range(nterm-1,0,-1):
            log1pxox = 1.e0 / j - x * log1pxox

    return log1pxox

class ModelMaker(object):
    def __init__(self, site_lat, site_lon, model_date, **kwargs):
        self.site_lat = site_lat
        self.site_lon = site_lon

        if(self.site_lon < 0.0):
            self.site_lon = self.site_lon + 360.0

        if(self.site_lon > 180.0):
            self.site_lon_180 = self.site_lon-360.0
        else:
            self.site_lon_180 = self.site_lon

        # If only a date is specified, just use noon time
        if len(model_date) == 10:
          self.model_dt = datetime.strptime(model_date, "%Y-%m-%d") + timedelta(hours=12)
        elif len(model_date) >= 19:
          self.model_dt = datetime.strptime(model_date, "%Y-%m-%dT%H:%M:%S")
        else:
          raise Exception("Unknown datetime format: %s" % model_date)

        # Load data from file
        self.ncep_data = self.load_ncep_data(**kwargs)

        # These are created as needed
        self.interp_data = None
        self.model = None

    def epoch_hours(self, ncep_obj):
        """Convert site date into number of hours at noon since epoch used by NCEP
            which is hours since 1-1-1 00:00:0.0 for NetCDF3 files and
            1800-01-01 00:00:0.0 for NetCDF4 files
        
            - self.site_lon_180/15.0 -- Account for offset of sun due to longitude 
            + 48 -- Epoch time is off by 2 days if you just convert
                    number of hours + 1,1,1 you can see this, I suppose due to being 
                    based on julian time? At least it was 2 days off for the first
                    Entry in the 2012 files..."""

        if(ncep_obj.file_format.find("NETCDF4") == 0):
            epoch_start = datetime(1800,1,1,0)
            correction = 0
        else:
            epoch_start = datetime(1,1,1,0)
            correction = 48
        epoch_hours = (self.model_dt - epoch_start).total_seconds()/60.0/60.0 \
            - self.site_lon_180/15.0 + correction

        return epoch_hours

    def load_ncep_data(self, air_temp_file=None, spec_hum_file=None, trop_pres_file=None, geo_height_file=None, surface_pres_file=None):
        """Loads NCEP Reanalysis data into NetCdf objects. Returns a named tuple with objects."""
        filenames = (air_temp_file, spec_hum_file, trop_pres_file, geo_height_file, surface_pres_file)
        dflt_search_globs = ("air.*.nc", "shum.*.nc", "pres.tropp.*.nc", "hgt.*.nc", "pres.sfc.*.nc")

        data_list = []
        for data_filename, srch_glob, file_type_name in zip(filenames, dflt_search_globs, NCEPData._fields):
            if not data_filename:
                srch_res = glob(srch_glob)
                if len(srch_res) != 1:
                    raise Exception("For file type: %s expected to find 1 NCEP file, instead found: %s, using glob: %s" % (file_type_name, srch_res, srch_glob))
                data_filename = srch_res[0]
            logger.debug("Loading NCEP file for %s: %s" % (file_type_name, data_filename))
    
            # Add filename for use elsewhere
            ncep_obj = NetCDFFile(data_filename, "r")
            ncep_obj.__dict__["filename"] = data_filename
    
            data_list.append(ncep_obj)

        return NCEPData(*data_list)

    def interp_to_site(self, ncep_obj, global_data_name):
        def frac_lat(lat_arr):
            return (self.site_lat - lat_arr[0]) / (lat_arr[1] - lat_arr[0])

        def frac_lon(lon_arr):
            return (self.site_lon - lon_arr[0]) / (lon_arr[1] - lon_arr[0])

        logger.debug("Loading data for %s interpolation from %s:" % (global_data_name, ncep_obj.filepath()))

        lat_file = ncep_obj.variables["lat"]
        logger.debug("  latitude: %s" % str(lat_file.shape))

        lon_file = ncep_obj.variables["lon"]
        logger.debug("  longitude: %s" % str(lon_file.shape))

        time_file = ncep_obj.variables["time"]
        logger.debug("  time: %s" % str(time_file.shape))

        global_data = ncep_obj.variables[global_data_name]
        logger.debug("  %s: %s" % (global_data_name, str(global_data.shape)))

        dlon = lon_file[-1] - lon_file[0]
        if (self.site_lat - lat_file[-1]) * (self.site_lat - lat_file[0]) > 0 or \
                (self.site_lon - lon_file[-1]) * (self.site_lon - lon_file[0]) * dlon > 0:
            raise Exeception("The data object for %s does not cover this site" % global_data_name)

        logger.debug("Interpolating %s data" % global_data_name)

        if hasattr(global_data, "scale_factor"):
            # Old NetCDF3 file with integers scaled and with offset
            global_data = global_data[:] * global_data.scale_factor + global_data.add_offset

        # Interpolate to lat/lon of site and to local noon
        frac_time = (self.epoch_hours(ncep_obj) - time_file[0]) / (time_file[1] - time_file[0])

        if len(global_data.shape) == 4:
            interp_loc = bilinear_interp2(global_data, frac_lon(lon_file), frac_lat(lat_file))
            interp_time = linear_interp1(interp_loc, frac_time)
        else:
            interp_loc = bilinear_interp1(global_data, frac_lon(lon_file), frac_lat(lat_file))
            interp_time = linear_interp0(interp_loc, frac_time)

        return interp_time

    def interpolate_ncep_to_site(self):
        logger.debug("Interpolating file data")

        # Check that the full number of levels have been downloaded
        # Geopotential & Air Temperature files
        if(self.ncep_data.air_temp.variables["level"].shape < 17):
            raise Exception("You must download all 17 levels of air temperature data from the NCEP/NCAR website.")

        if(self.ncep_data.geo_height.variables["level"].shape < 17):
            raise Exception("You must download all 17 levels of geopotential height data from the NCEP/NCAR website.")

        if(self.ncep_data.spec_hum.variables["level"].shape < 8):
            raise Exception("You must download all 8 levels of specific humidity data from the NCEP/NCAR website.")

        # Air Temperature interpolated to Lat/Long of site:
        air_temp_data = self.interp_to_site(self.ncep_data.air_temp, "air")

        # Geopotential Height interpolated to Lat/Long of site:
        geo_height_data = self.interp_to_site(self.ncep_data.geo_height, "hgt")

        # Convert m to km
        geo_height_data /= 1000.0

        # Tropopause Pressure interpolated to Lat/Long of site:
        trop_pres_data = self.interp_to_site(self.ncep_data.trop_pres, "pres")

        # Specific Humidity interpolated to Lat/Long of site:
        # Convert the spec hum data into vmr units
        spec_hum_data = self.interp_to_site(self.ncep_data.spec_hum, "shum")
        spec_hum_data *= 28.964/18.02

        # Surface pressure interpolated to Lat/Lon of site
        surf_pres_data = self.interp_to_site(self.ncep_data.surface_pres, "pres")

        self.interp_data = NCEPData(air_temp_data, spec_hum_data, trop_pres_data, geo_height_data, surf_pres_data)
    
        return self.interp_data
    
    def create_model(self):
        # Creates a model data structure
        #     Site_Lat    ; The latitude of the site
        #     Lev_AT      ; The pressure levels on which the data are tabulated
        #     sat         ; Site Noon Atmospheric Temperature profile (vector)
        #     sgh         ; Site Noon Geometric Height profile in km (vector)
        #     ssh         ; Site Noon Specific Humidity profile (vector)
        #     stp         ; Site Noon Tropopause Pressure (scalar)
        
        # Calculate interpolated data if it has not been loaded yet
        if not self.interp_data:
            self.interpolate_ncep_to_site()

        lev_at = self.ncep_data.air_temp.variables["level"]
        sat = self.interp_data.air_temp
        sgh = self.interp_data.geo_height
        ssh = self.interp_data.spec_hum
        stp = self.interp_data.trop_pres

        logger.debug("Creating model data")
        p_out = []
        temp_out = []
        height_out = []
        h2o_out = []
        surf_pres_out = self.interp_data.surface_pres

        # Export the Pressure, Temp and SHum for lower levels (1000 to 300 mbar)
        for k in range(0, len(ssh)):
            if (ssh[k] < 4.0e-06 and k > 0):
                logger.debug("Replacing insufficient SHum %f %f" (lev_at[k], ssh[k]))
                ssh[k] = sqrt(4.0e-06 * ssh[k-1])

            p_out.append(lev_at[k])            
            temp_out.append(sat[k])
            height_out.append(sgh[k])
            h2o_out.append(ssh[k])

        # Export Pressure and Temp for middle levels (250 to 10 mbar)
        # which have no SHum reanalysis.
        dsh = ssh[len(ssh)-1]
        lstp=log10(stp)
        for k in range(len(ssh), lev_at.shape[0]):
            zz = log10(lev_at[k])  # log10[pressure]
            strat_h2o=7.5E-06*exp(-3.0*(zz/lstp)**2)
            dsh = dsh * pow(lev_at[k] / lev_at[k-1], 5.5 - lev_at[k] / 100)
            svp = svp_wv_over_ice(sat[k])
            satvmr = svp / lev_at[k]   # vmr of saturated WV at T/P

            if ( dsh > satvmr ):
                dsh=satvmr
            if (lev_at[k] < stp / 100.0):
                dsh = sqrt(dsh * strat_h2o)  # above the trop

            p_out.append(lev_at[k])
            temp_out.append(sat[k])
            height_out.append(sgh[k])
            h2o_out.append(max([dsh,strat_h2o]))

        # Get the difference between the USSA and given site temperature at 10 mbar,
        delta_t = sat[16] - t_ussa[0]

        # Export the P-T profile above 10mbar
        for k in range(1, 8):
            delta_t = delta_t/2
            zz = log10(p_ussa[k])  # log10[pressure]
            strat_h2o=7.5E-06*exp(-0.25*zz**2)

            p_out.append(p_ussa[k])
            temp_out.append(t_ussa[k]+delta_t)
            height_out.append(z_ussa[k])
            h2o_out.append(strat_h2o)

        units = ("Pa", "Kelvin", "km", "vmr", "Pa")
        data_name = ("Pressure", "Temperature", "Height", "H2O", "Surface_Pressure")
        data_in = (p_out, temp_out, height_out, h2o_out, surf_pres_out)

        data_list = []
        for u, nm, dat in zip(units, data_name, data_in):
            data_list.append(ModelData(nm, u, array(dat)))
        
        self.model = Model(*data_list)
        return self.model

    def interpolate_model_to_grid(self, z_arr=Z_GRID, w_arr=MEAN_MOL_WGT):
        """Interpolate model values to a specified (or default) height grid (in km) with associated mean
        molecular weight values. 

        This routine is based off of readmodFC.f from gsetup
        """

        if not self.model:
            self.create_model()

        nlev = z_arr.shape[0]
        ninlvl = self.model.pressure.data.shape[0]

        temp_interp = zeros(nlev, dtype=float)
        press_interp = zeros(nlev, dtype=float)
        dens_interp = zeros(nlev, dtype=float)
        h2o_interp = zeros(nlev, dtype=float)

        # Constants written into the .mod file by modmaker
        # Radius in modmaker is actually inconsistent in that they
        # use two different values with the value that would have
        # been written to the file being 6378.00 instead of the
        # constant used below
        radius = EARTH_RADIUS
        ztrop_gct=0.0  # initial values
        ztrop_ncep=0.0 # initial values
        ecc2 = 6.0e-5
        tlat = self.site_lat
        gs = 9.81
        zold = self.interp_data.geo_height[0]
        pfact = 1013.25
        ptrop_ncep = self.interp_data.trop_pres/10**2  # mbar

        #print radius,ecc2, tlat,gs,zold,pfact,ptrop_ncep
        hold=zold/(1+zold/radius)  # Convert geometric to geopotential altitude
        roc=radius                 # don't correct for Earth's ecentricity
                                   # OK for ground based & mid latitude

        inpress = self.model.pressure.data
        intemp = self.model.temperature.data
        inheight = self.model.height.data
        inh2ovmr = self.model.h2o_vmr.data

        pold=inpress[0]
        told=intemp[0]
        zold=inheight[0]
        h2oold=inh2ovmr[0]
        #print "nlev=",nlev,pold,told,zold

        if (nlev > 1):
            pnew=inpress[1]
            tnew=intemp[1]
            znew=inheight[1]
            h2onew=inh2ovmr[1]
            lr=-5.0
            zlr=0.0
            ii=2
        else:
            hnew=1.E+36
            pnew=pold
            tnew=told
            znew=zold
            h2onew=h2oold

        #print hold,pold,told
        x=tnew/told-1
        hnew=hold+U_GAS_CONST*told*log(pold/pnew)/gs/w_arr[0]/log1pxox(x)
        #print 'hnew=',hnew,pold,pnew,told,tnew,x
        #print 'hnew=',hnew,hnew/(1-hnew/radius),pnew,tnew

        iztrop=0.0
        for k in range(nlev):       # loop over levels
            hk=z_arr[k]/(1+z_arr[k]/radius)  # Convert from geometric to geopotential
            while (hk > hnew):      # hold & hnew are both below hk; read another record
                pold=pnew
                told=tnew
                zold=znew
                hold=hnew
                lrwas=lr
                zlrwas=zlr
                h2oold=h2onew 
                ii=ii+1
                if(ii > ninlvl):
                    logger.debug('ii, ninlvl, modname z(k) zold znew = %f, %d, %d, %f, %f, %f' % \
                        (ii,ninlvl,ninlvl,z(k),hold/(1-hold/radius),hnew/(1-hnew/radius)))
                    raise Exception('Warning: Model levels do not extend high enough Model levels do not extend high enough')
                pnew=inpress[ii]
                tnew=intemp[ii]
                znew=inheight[ii]

                # Compute tropopause altitude
                lr=(tnew-told)/(znew-zold) # Lapse Rate
                zlr=0.5*(zold+znew)        # Altitude at which lapse rate = lr

                # Find first instance of lapse-rate exceeding -2K/km
                # Only valid for Earth:
                #print k, zold, told, zlr, lr
                if(iztrop == 0.0 and zold > 5 and lr > -2.0):
                    iztrop=1.0 # do not compute the trop altitude more than once
                    ztrop_gct=zlrwas+(zlr-zlrwas)*(-2-lrwas)/(lr-lrwas)
                    ztrop_gct=ztrop_gct/(1-ztrop_gct/radius)  # convert H to Z

                if(tnew < 0.0):
                    logger.debug("pnew, tnew = %f, %f" % pnew, tnew)
                    raise Exception('Model temperatures must be in Kelvin')

                x=tnew/told-1
                hnew=hold+U_GAS_CONST*told*log(pold/pnew)/gs/w_arr[k]/log1pxox(x)
                #print hnew,hnew/(1-hnew/radius),pnew,tnew
                h2onew=inh2ovmr[ii]

            #   z(k) is bracketed by zold & znew; proceed with interpolation
            h2ox=0.0
            if(hnew == hold):
                x = 0.0
            else:
                x=(tnew/told-1)*(hk-hold)/(hnew-hold)
                if(h2oold > 0.0):
                    h2ox=(h2onew/h2oold-1)*(hk-hold)/(hnew-hold)
            temp_interp[k] = told*(1+x)
            press_interp[k] = pold/pfact*exp((hold-hk)*w_arr[k]*gs*log1pxox(x)/U_GAS_CONST/told)
            #print 'p(k)=', press_interp[k],told,tnew,hk,hold,hnew,U_GAS_CONST,w_arr[k],gs,x
            dens_interp[k] = 0.101325*AVAGADRO*press_interp[k]/temp_interp[k]/U_GAS_CONST  #units of cm-3
            h2o_interp[k] = h2oold*(1+h2ox) #RAW Linear interpolation of H2O

        #print zold,pold,told,z_arr[k],temp_interp[k],press_interp[k],h2o_interp[k]

        # Convert NCEP tropopause pressure to altitude
        if(ptrop_ncep > 0 and  nlev > 1):
            for k in range(2,nlev): # loop over levels
                if(press_interp[k] < ptrop_ncep / pfact):
                  break

            fr=log(ptrop_ncep/pfact/press_interp[k])/log(press_interp[k-1]/press_interp[k])
            ztrop_ncep=fr*z_arr[k]+(1-fr)*z_arr[k]
        else:
            ztrop_ncep=0.0

        # Convert back to milibars
        press_interp *= pfact

        # Replace model values
        data_models = (self.model.pressure, self.model.temperature, self.model.height, self.model.h2o_vmr, self.model.surface_pressure)
        new_data = (press_interp, temp_interp, z_arr, h2o_interp, self.model.surface_pressure.data)

        data_list = []
        for mod, data in zip(data_models, new_data):
            data_list.append(mod._replace(data=data))

        self.model = Model(*data_list)
        return self.model

    def write_oco_ascii(self, filename=None):
        if not self.model:
            self.create_model()

        fn_date_str = self.model_dt.strftime("%Y%m%d%H%M")
        id_date_str = self.model_dt.strftime("%Y-%m-%dT%H:%M:%S")
        if not filename:
            filename = "model_%s.dat" % (fn_date_str)

        nlev = self.model.pressure.data.shape[0]

        # Convert pressure from mbar to Pascals
        press_mbar = self.model.pressure 
        press_pa = press_mbar._replace(data=press_mbar.data * 100.0)
        data_in_order = (press_pa, self.model.height, self.model.temperature, self.model.h2o_vmr)
        data_arr = zeros((nlev, len(data_in_order)), dtype=float)

        for idx in range(len(data_in_order)):
            # Reverse data to be increasing pressure order
            data_arr[:, idx] = data_in_order[idx].data[::-1]

        # Put import here to not make this code rely on this module
        from full_physics.oco_matrix import OcoMatrix

        out_mat = OcoMatrix()

        out_mat.file_id = "Atmosphere Model Created from NCEP data interpolated to latitude: %f, longitude: %f on %s" % (self.site_lat, self.site_lon, id_date_str)

        out_mat.units = [] 
        out_mat.labels = []
        for val in data_in_order:
            out_mat.units.append(val.units)
            if val.name == "Temperature":
                out_mat.labels.append("T")
            else:
                out_mat.labels.append(val.name)

        out_mat.header["Surface_Pressure"] = self.model.surface_pressure.data
        out_mat.header["Surface_Pressure_Units"] = self.model.surface_pressure.units

        out_mat.data = data_arr

        logger.debug("Writing to ASCII file: %s" % filename)
        out_mat.write(filename)

def standalone_main():
    parser = OptionParser(usage="usage: %prog --lat <site_latitude> --lon <site_longitude> -d <model_date>")

    parser.add_option( "--lat", "--latitude", dest="latitude",
                       metavar="NUMBER", type="float",
                       help="latitude of location where data is interpolated")

    parser.add_option( "--lon", "--longitude", dest="longitude",
                       metavar="NUMBER", type="float",
                       help="longitude of location where data is interpolated")

    parser.add_option( "-d", "--model_date", dest="model_date",
                       metavar="DATETIME",
                       help="date string representing the date of model in the format YYYY-MM-DD")

    parser.add_option( "-i", "--interp", dest="interp",
                       action="store_true",
                       default=False,
                       help="interpolate to a higher number of levels grid")

    parser.add_option( "-o", "--output_file", dest="output_file",
                       metavar="FILE",
                       help="filename where data is output")

    parser.add_option( "--air_temp", dest="air_temp_file",
                       metavar="FILE",
                       help="NCEP air temperature file")

    parser.add_option( "--spec_hum", dest="spec_hum_file",
                       metavar="FILE",
                       help="NCEP specific humidity file")

    parser.add_option( "--trop_pres", dest="trop_pres_file",
                       metavar="FILE",
                       help="NCEP troposphere pressure file")

    parser.add_option( "--geo_height", dest="geo_height_file",
                       metavar="FILE",
                       help="NCEP geopotential height file")

    parser.add_option( "--surface_pres", dest="surface_pres_file",
                       metavar="FILE",
                       help="NCEP surface pressure file")

    parser.add_option( "-v", "--verbose", dest="verbose",
                       action="store_true",
                       default=False,
                       help="enable verbose informational reporting")


    # Parse command line arguments
    (options, args) = parser.parse_args()

    if not options.latitude or not options.longitude:
        parser.error("Must specify latitude and longitude of interpolation location")

    if not options.model_date:
        parser.error("Must specify model date string")

    # Set up logging
    if options.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    console = logging.StreamHandler()
    logger.addHandler(console)

    # Load NCEP reanalysis files, optionally specified on the command line
    ncep_data_files = {}
    for opt_nm in ("air_temp_file", "spec_hum_file", "trop_pres_file", "geo_height_file", "surface_pres_file"):
        ncep_data_files[opt_nm] = getattr(options, opt_nm)
    mm = ModelMaker(options.latitude, options.longitude, options.model_date, **ncep_data_files)

    if options.interp:
        mm.interpolate_model_to_grid()
    
    mm.write_oco_ascii(options.output_file)


if __name__ == "__main__":
    standalone_main()
