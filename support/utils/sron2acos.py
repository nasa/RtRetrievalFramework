#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import sys
import h5py
from numpy import *
from optparse import OptionParser

def parseSpectraFile(spec_file, f, name, ll, offset):
    f_band = open(spec_file,'r')
    band = f_band.readlines()
    numberOfSoundings = len(band)-5 # Usually, first 5 lines are comments as to the file structure, be careful here! 
    line10 = band[5].split()       # Quick& dirty, split line 10 to get real spectra size
    specLength = old_div((len(line10)-1),4)  # One line contains filename + two pairs of wavenumbers and spectra, i.e. (length-1) divided by 4 gives length of a single spectrum
    print('Parsing spectra, might take some time...')
    spec = f.create_dataset(name, (numberOfSoundings,2,ll), 'f')
    for i in range(5,len(band)):
        line_band = band[i].split()
        # set P polarization spectra
        startIndex = 1+2*specLength
        for ii in range(startIndex, startIndex+specLength):
            spec[i-5, 0,ii-startIndex+offset] = float(line_band[ii])
        # set S polarization spectra
        startIndex = 1+3*specLength
        for ii in range(startIndex, startIndex+specLength):
            spec[i-5, 1,ii-startIndex+offset] = float(line_band[ii])
    f_band.close()
    return specLength, numberOfSoundings

def parseGeolocations(geoFile, f):
    geo = open(geoFile,'r')
    geo_line = geo.readlines()
    numberOfSoundings = len(geo_line)-9 # Usually, first 9 lines are comments as to the file structure, be careful here! 
    print('Creating Footprint Geometries for ', numberOfSoundings, ' soundings')
    # Geometries and Gain settings...
    lat = f.create_dataset('FootprintGeometry/footprint_latitude', (numberOfSoundings,3,2), 'f')
    lon = f.create_dataset('FootprintGeometry/footprint_longitude', (numberOfSoundings,3,2), 'f')
    alt = f.create_dataset('FootprintGeometry/footprint_altitude', (numberOfSoundings,3,2), 'f')
    sza = f.create_dataset('FootprintGeometry/footprint_solar_zenith', (numberOfSoundings,3,2), 'f')
    saa = f.create_dataset('FootprintGeometry/footprint_solar_azimuth', (numberOfSoundings,3,2), 'f')
    lza = f.create_dataset('FootprintGeometry/footprint_zenith', (numberOfSoundings,3,2), 'f')
    laa = f.create_dataset('FootprintGeometry/footprint_azimuth', (numberOfSoundings,3,2), 'f')
    time = f.create_dataset('FootprintGeometry/footprint_time', (numberOfSoundings,3,2), 'd')
    gain = f.create_dataset('SoundingHeader/gain_swir',  (numberOfSoundings,2), '|S6')
    # Determine noise values...
    cnv1 = ones((len(geo_line)-9,2,1805), dtype='f' )
    cnv2 = ones((len(geo_line)-9,2,3508), dtype='f' )
    cnv3 = ones((len(geo_line)-9,2,2005), dtype='f' )
    noise_1 = 1.5e-9*ones((len(geo_line)-9,2), dtype='f' )
    noise_2 = 1.5e-9*ones((len(geo_line)-9,2), dtype='f' )
    noise_3 = 1.5e-9*ones((len(geo_line)-9,2), dtype='f' )
    f['InstrumentHeader/cnv_coef_highgain_o2'] =  cnv1
    f['InstrumentHeader/cnv_coef_medgain_o2'] =  cnv1
    f['InstrumentHeader/cnv_coef_highgain_weak_co2'] =  cnv2
    f['InstrumentHeader/cnv_coef_medgain_weak_co2'] =  cnv2
    f['InstrumentHeader/cnv_coef_highgain_strong_co2'] =  cnv3
    f['InstrumentHeader/cnv_coef_medgain_strong_co2'] =  cnv3
    f['SoundingSpectra/noise_o2'] = noise_1
    f['SoundingSpectra/noise_weak_co2'] = noise_2
    f['SoundingSpectra/noise_strong_co2'] = noise_3
    
    wn_coeff = f.create_dataset('SoundingHeader/wavenumber_coefficients',  (numberOfSoundings,3,2,2), 'd')
    velocity = zeros( (len(geo_line)-9), dtype='f' )
    f['SpacecraftGeometry/relative_velocity'] =  velocity
    stokes = f.create_dataset('FootprintGeometry/footprint_stokes_coefficients', (numberOfSoundings,3,2,4), 'f')
    psurf = zeros((numberOfSoundings,1))
    sid = ones( (len(geo_line)-9), dtype=int64 )
    wn_coeff[:,0,:,0] = 12870.084
    wn_coeff[:,1,:,0] = 5750.183
    wn_coeff[:,2,:,0] = 4750.125
    wn_coeff[:,:,:,1] = 0.19949289000000001
    for i in range(9,len(geo_line)):
        geo_line_split = geo_line[i].split()
        splitter = geo_line_split[0].split('_')
        sounding_id = float('2007'+splitter[7]+splitter[8]+'0000'+str(i-9).zfill(2))
     #   print sounding_id
        sid[i-9] = sounding_id
        lat[i-9,:,:] = float(geo_line_split[3]);
        lon[i-9,:,:] = float(geo_line_split[4]);
        alt[i-9,:,:] = float(geo_line_split[1]);
        sza[i-9,:,:] = float(geo_line_split[5]);
        saa[i-9,:,:] = float(geo_line_split[6]);
        lza[i-9,:,:] = float(geo_line_split[7]);
        laa[i-9,:,:] = float(geo_line_split[8]);
        gain[i-9,:]='H     '
        time[i-9,:,:] = 4.323879e8 #fill dummy value for time here, don't know whether it is really used...
        stokes[i-9,:,:,:] = [1,0,0,0]
        psurf[i-9] = float(geo_line_split[2]);
    f['SoundingHeader/sounding_id'] = sid
    f['SoundingGeometry/sounding_id'] = sid
    f['FootprintGeometry/footprint_time_tai93'] = time
    geo.close()
    return numberOfSoundings, psurf

def parseMetFile(metFile, f, psurf):
    met = open(metFile,'r')
    met_line = met.readlines()
    print(len(met_line))
    numberOfSoundings = len(met_line)-4 # Usually, first 4 lines are comments as to the file structure, be careful here! 
    fieldLength = old_div((len(met_line[4].split())-1),3) # Field length
    print('Reading Met files...')
    print('Number of soundings in met file: ', numberOfSoundings)
    print('Number of atmospheric layers: ', fieldLength)
    press_level1 = f.create_dataset('ecmwf/specific_humidity_pressures', (numberOfSoundings,3,2, fieldLength), 'f')
    press_level2 = f.create_dataset('ecmwf/temperature_pressures', (numberOfSoundings,3,2, fieldLength), 'f')
    surfacePressure = f.create_dataset('ecmwf/surface_pressure', (numberOfSoundings,3,2), 'f')
    temperature =  f.create_dataset('ecmwf/temperature', (numberOfSoundings,3,2, fieldLength), 'f')
    specHumidity = f.create_dataset('ecmwf/specific_humidity', (numberOfSoundings,3,2, fieldLength), 'f')
    wind_u = f.create_dataset('ecmwf/windspeed_u', (numberOfSoundings,3,2), 'f')
    wind_v = f.create_dataset('ecmwf/windspeed_v', (numberOfSoundings,3,2), 'f')
    wind_u[:,:,:]=0
    wind_v[:,:,:]=0
    for i in range(4,len(met_line)):
        surfacePressure[i-4,:,:] = 100*psurf[i-4]
        met_line_split = met_line[i].split()
        startIndex = 1+0*fieldLength
        # Writing pressure grids in Pa
        for ii in range(startIndex, startIndex+fieldLength):
            pr = 100*float(met_line_split[ii])
         #   print pr
            if pr<100:
                pr=100
            press_level1[i-4, :,:,ii-startIndex] = pr
            press_level2[i-4, :,:,ii-startIndex] = pr
        startIndex = 1+1*fieldLength
        # Writing the temperature array
        for ii in range(startIndex, startIndex+fieldLength):
            temperature[i-4, :,:,ii-startIndex] = float(met_line_split[ii])
        startIndex = 1+2*fieldLength
        # read H2O dry air mole fraction and convert to specific humidity
        for ii in range(startIndex, startIndex+fieldLength):
            h2o_vmr = float(met_line_split[ii])
            specHumidity[i-4, :,:,ii-startIndex] =0.622*h2o_vmr/(1+h2o_vmr)
            
    met.close()
    
def standalone_main():
    parser = OptionParser(usage="usage: %prog --geo geoLocationsFile.dat --met metFieldsFile.dat --spec_1 o2A_band_spectra_file --spec_2 wCO2_band_spectra_file --spec_3 sCO2_band_spectra_file -l l1B_outputFile -e ECMWF_outputFile ")

    parser.add_option( "--geo", "--geolocation_file", dest="geo_file",
                       metavar="FILE",
                       help="text file with geolocation specification from the SRON orbit simulator")
    parser.add_option( "--met", "--meteorology_file", dest="met_file",
                       metavar="FILE",
                       help="text file with meteorology profiles used in the SRON orbit simulator")
    parser.add_option( "--spec_1", "--o2_A_band_spectraFile", dest="spec_band1",
                       metavar="FILE",
                       help="text file with O2A band (0.76 micron) spectra generated by the SRON orbit simulator")
    parser.add_option( "--spec_2", "--weak_co2_band_spectraFile", dest="spec_band2",
                       metavar="FILE",
                       help="text file with weak CO2 band (1.6 micron) spectra generated by the SRON orbit simulator")
    parser.add_option( "--spec_3", "--strong_co2_band_spectraFile", dest="spec_band3",
                       metavar="FILE",
                       help="text file strong CO2 band (2 micron) spectra generated by the SRON orbit simulator")
    parser.add_option( "-l", "--l1b_output_file", dest="l1b_output_file",
                       metavar="FILE",
                       default="l1b_fromSRON.h5",
                       help="generated l1b HDF5 output filename (according to ACOS fileformat definitions)")
    parser.add_option( "-e", "--ecm_output_file", dest="ecm_output_file",
                       metavar="FILE",
                       default="ecm_fromSRON.h5",
                       help="generated ECMWF HDF5 output filename (according to ACOS fileformat definitions)")


    # Parse command line arguments
    (options, args) = parser.parse_args()

    if options.geo_file == None:
        parser.error('Geolocation file must be specified')

    if options.met_file == None:
        parser.error('Meterology input file must be specified')
    
    # Create L1B file
    f = h5py.File(options.l1b_output_file,'w')
    # create necessary groups first:
    f.create_group("SoundingSpectra")
    f.create_group("SoundingGeometry")
    f.create_group("FootprintGeometry")
    f.create_group("SoundingHeader")
    f.create_group("SpacecraftGeometry")
    f.create_group("InstrumentHeader")
    
    # Create ECMWF output file
    f_m =   h5py.File(options.ecm_output_file,'w')
    f_m.create_group("ecmwf")

    # Read gelocations and write them into HDF5 file:
    numberOfSoundings, psurf = parseGeolocations(options.geo_file, f)
    # Read meteorology files and write them into HDF5 file:
    parseMetFile(options.met_file, f_m, psurf)
    
    # parse all spectra files into HDF5 file:
    specLength1, numberOfSoundings1 =  parseSpectraFile(options.spec_band1,f,  '/SoundingSpectra/radiance_o2', 1805, 401)
    specLength2, numberOfSoundings2 =  parseSpectraFile(options.spec_band2,f,  '/SoundingSpectra/radiance_weak_co2', 3508, 2075)
    specLength3, numberOfSoundings3 =  parseSpectraFile(options.spec_band3,f,  '/SoundingSpectra/radiance_strong_co2', 2005, 281)
    
    if numberOfSoundings2 == numberOfSoundings and numberOfSoundings3 == numberOfSoundings and numberOfSoundings1 == numberOfSoundings:
        print("Number of soundings consistent, namely: ", numberOfSoundings)
        print("O2A band spectra length: " , specLength1)
        print("wCO2 band spectra length: " , specLength2)
        print("sCO2 band spectra length: " , specLength3)
    else:
        print("Error: Number of soundings inconsistent, namely: \n")
        print("Geolocations: ",  numberOfSoundings)
        print("O2A band file: " , numberOfSoundings1)
        print("wCO2 band file: " , numberOfSoundings2)
        print("sCO2 band file: " , numberOfSoundings3)
        sys.exit(2)
    
   
    f.close()
    f_m.close()
    
if __name__ == "__main__":
    standalone_main()
