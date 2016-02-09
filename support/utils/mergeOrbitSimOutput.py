#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
import h5py
import sys
from numpy import *

import glob
import full_physics.acos_file as acos_file
from optparse import OptionParser

def mergeFiles(options):
    l1b_file = options.l1b_file 
    print('Using '+  l1b_file + ' to extract header information')
    l1b_obj = acos_file.L1B(l1b_file)
    l2_obj = h5py.File(options.l2_file, "r+")
    try:
     #   print 'trying hard...'
        l2_obj.create_group("SoundingHeader")
        l2_obj.create_group("SoundingGeometry")
    except:
        print("Sounding header already exists in the L2 file!")
        sys.exit(0)
    exposure_index = l2_obj["RetrievalResults/exposure_index"][:]-1
    sounding_id_reference = l2_obj["/RetrievalResults/sounding_id_reference"][:]
    time_string = l1b_obj["SoundingHeader/exposure_start_time_string"][:][exposure_index]
    try:    
        start_time = l1b_obj["SoundingHeader/exposure_start_time_tai93"][:][exposure_index]
    except:
        print("exposure_start_time not found (probably simulator file?), skipping this...")
   
    # write few entries into L2 files...
    l2_obj["SoundingHeader/exposure_start_time_string"] = time_string
    try:    
        l2_obj["SoundingHeader/exposure_start_time_tai93"] = start_time
    except:
        print("time missing")
    l2_obj["SoundingHeader/sounding_id"] = l1b_obj["SoundingHeader/sounding_id"][:][exposure_index]
    l2_obj["SoundingGeometry/sounding_longitude"] = l1b_obj["SoundingGeometry/sounding_longitude"][:][exposure_index]
    l2_obj["SoundingGeometry/sounding_latitude"] = l1b_obj["SoundingGeometry/sounding_latitude"][:][exposure_index]
    l2_obj["SoundingGeometry/sounding_solar_zenith"] =  l1b_obj["SoundingGeometry/sounding_solar_zenith"][:][exposure_index]
    l2_obj["SoundingGeometry/sounding_zenith"] =  l1b_obj["SoundingGeometry/sounding_zenith"][:][exposure_index]
    l2_obj["SoundingGeometry/sounding_solar_azimuth"] = l1b_obj["SoundingGeometry/sounding_solar_azimuth"][:][exposure_index]
    l2_obj["SoundingGeometry/sounding_azimuth"] =l1b_obj["SoundingGeometry/sounding_azimuth"][:][exposure_index]
    l2_obj["SoundingGeometry/sounding_altitude"] =l1b_obj["FootprintGeometry/footprint_altitude"][:][exposure_index]
    # sanity check:
    # print sounding_id_reference[0], sounding_id[0]
    
    ###################################################################################
    # Proceed with the log file:
    l2_obj.create_group("truth")
    l2_obj.create_group("truth/log")
    l2_obj.create_group("truth/scene")
    l2_obj.create_group("truth/simulation")
    logData = loadtxt(options.log_file, comments='#')
    l2_obj["truth/log/frame"] = int32(logData[exposure_index,0])
    l2_obj["truth/log/sounding_id"] = int64(logData[exposure_index,1])
    l2_obj["truth/log/surface_model"] = logData[exposure_index,2]
    l2_obj["truth/log/albedo_1A"] = logData[exposure_index,3]
    l2_obj["truth/log/albedo_1B"] = logData[exposure_index,4]
    l2_obj["truth/log/albedo_2A"] = logData[exposure_index,5]
    l2_obj["truth/log/albedo_2B"] = logData[exposure_index,6]
    l2_obj["truth/log/albedo_3A"] = logData[exposure_index,7]
    l2_obj["truth/log/albedo_3B"] = logData[exposure_index,8]
    l2_obj["truth/log/tau_water_1"] = logData[exposure_index,9]
    l2_obj["truth/log/tau_ice_1"] = logData[exposure_index,10]
    l2_obj["truth/log/tau_aerosol_1"] = logData[exposure_index,11]
    l2_obj["truth/log/tau_water_2"] = logData[exposure_index,12]
    l2_obj["truth/log/tau_ice_2"] = logData[exposure_index,13]
    l2_obj["truth/log/tau_aerosol_2"] = logData[exposure_index,14]
    l2_obj["truth/log/tau_water_3"] = logData[exposure_index,15]
    l2_obj["truth/log/tau_ice_3"] = logData[exposure_index,16]
    l2_obj["truth/log/tau_aerosol_3"] = logData[exposure_index,17]
    l2_obj["truth/log/logFile"] = options.log_file
    ###################################################################################
    # Proceed with the scene log file:
    scene = loadtxt(options.scene_file, comments='#')
    l2_obj["truth/scene/sceneFile"] =options.scene_file
    l2_obj["truth/scene/frame"] = int32(scene[exposure_index,0])
    l2_obj["truth/scene/sounding_id"] = int64(scene[exposure_index,1])
    l2_obj["truth/scene/sounding"] = int32(scene[exposure_index,2])
    l2_obj["truth/scene/profile_index"] = scene[exposure_index,3]
    l2_obj["truth/scene/sat_latitude"] = scene[exposure_index,4]
    l2_obj["truth/scene/sat_longitude"] = scene[exposure_index,5]
    l2_obj["truth/scene/fov_latitude"] = scene[exposure_index,6]
    l2_obj["truth/scene/fov_longitude"] = scene[exposure_index,7]
    l2_obj["truth/scene/sun_latitude"] = scene[exposure_index,8]
    l2_obj["truth/scene/sun_longitude"] = scene[exposure_index,9]
    l2_obj["truth/scene/XCO2"] = scene[exposure_index,10]
    l2_obj["truth/scene/psurf"] = scene[exposure_index,11]
    l2_obj["truth/scene/mean_temperature"] = scene[exposure_index,12]
    l2_obj["truth/scene/cloud_water"] = scene[exposure_index,13]   
    l2_obj["truth/scene/cloud_ice"] = scene[exposure_index,14]
    l2_obj["truth/scene/scattering_height_water"] = scene[exposure_index,15]
    l2_obj["truth/scene/scattering_height_ice"] = scene[exposure_index,16]
    l2_obj["truth/scene/aerosol_1"] = scene[exposure_index,17]
    l2_obj["truth/scene/aerosol_2"] = scene[exposure_index,18]
    l2_obj["truth/scene/aerosol_3"] = scene[exposure_index,19]
    l2_obj["truth/scene/aerosol_4"] = scene[exposure_index,20]
    l2_obj["truth/scene/aerosol_5"] = scene[exposure_index,21]
    l2_obj["truth/scene/aerosol_6"] = scene[exposure_index,22]
    l2_obj["truth/scene/aerosol_7"] = scene[exposure_index,23]
    l2_obj["truth/scene/sun_zenith"] = scene[exposure_index,24]
    l2_obj["truth/scene/sat_zenith"] = scene[exposure_index,25]
    l2_obj["truth/scene/IGBP_index"] = scene[exposure_index,26]
    l2_obj["truth/scene/surface_windSpeed"] = scene[exposure_index,27]
    l2_obj["truth/scene/moist_air"] = scene[exposure_index,28]
    l2_obj["truth/scene/h2o"] = scene[exposure_index,29]
    l2_obj["truth/scene/h2o_"] = scene[exposure_index,30]
# l2_obj["truth/scene/h2o_"].attrs['Units']='kg/m^2'
    ###################################################################################
    # Proceed with the detailed HDF scene file
    sim = h5py.File(options.logHDF_file, "r")
    l2_obj.create_group("truth/simulation/Aerosol")
    l2_obj.create_group("truth/simulation/Gas")
    l2_obj.create_group("truth/simulation/Surface")
    l2_obj.create_group("truth/simulation/Thermodynamic")
    l2_obj["truth/simulation/Aerosol/num_species"]=sim["Simulation/Aerosol/num_species"][:][exposure_index]
    l2_obj["truth/simulation/Aerosol/species_density"]=sim["Simulation/Aerosol/species_density"][:][exposure_index,:,:]
    l2_obj["truth/simulation/Aerosol/species_id"]=sim["Simulation/Aerosol/species_id"][:][exposure_index,:]
    l2_obj["truth/simulation/Gas/num_species"]=sim["Simulation/Gas/num_species"][:][exposure_index]
    l2_obj["truth/simulation/Gas/species_density"]=sim["Simulation/Gas/species_density"][:][exposure_index,:,:]
    l2_obj["truth/simulation/Gas/species_id"]=sim["Simulation/Gas/species_id"][:][exposure_index,:]

    l2_obj["truth/simulation/Surface/land_fraction"]=sim["Simulation/Surface/land_fraction"][:][exposure_index]
    l2_obj["truth/simulation/Thermodynamic/altitude_level"]= sim["Simulation/Thermodynamic/altitude_level"][:][exposure_index,:]
    l2_obj["truth/simulation/Thermodynamic/num_layers"]= sim["Simulation/Thermodynamic/num_layers"][:][exposure_index]
    l2_obj["truth/simulation/Thermodynamic/pressure_level"]= sim["Simulation/Thermodynamic/pressure_level"][:][exposure_index,:]
    l2_obj["truth/simulation/Thermodynamic/temperature_level"]= sim["Simulation/Thermodynamic/temperature_level"][:][exposure_index,:]
    ###################################################################################
    # Try to implement the averaging kernel correction (all done on the retrieval grid, not the simulator grid, to be though over)
    l2_obj.create_group("truth/co2_profile")
    press_level = l2_obj["RetrievalResults/vector_pressure_levels"][:]
    n_level = l2_obj["RetrievalResults/num_active_levels"][:]
    l2_obj["truth/co2_profile/p"] = l2_obj["RetrievalResults/vector_pressure_levels"]
    l2_obj["truth/co2_profile/n_levels"] = l2_obj["RetrievalResults/num_active_levels"]
    l2_obj["truth/co2_profile/co2_retrieved"] = l2_obj["RetrievalResults/co2_profile"]
    l2_obj["truth/co2_profile/co2_apriori"] = l2_obj["RetrievalResults/co2_profile_apriori"]
    # generate true profiles interpolated on the retrieval grid (has to be extrapolated at times, stupid!)
    co2_prior = l2_obj["RetrievalResults/co2_profile_apriori"]
    p_coarse = l2_obj["RetrievalResults/vector_pressure_levels"]
    n_coarse = l2_obj["RetrievalResults/num_active_levels"]
    n_fine = sim["Simulation/Thermodynamic/num_layers"][:][exposure_index]
    p = sim["Simulation/Thermodynamic/pressure_level"][:][exposure_index,:]
    gases = sim["Simulation/Gas/species_density"][:][exposure_index,:,:]
    AK_all = l2_obj["/RetrievalResults/averaging_kernel_matrix"]
    h =l2_obj["RetrievalResults/xco2_pressure_weighting_function"];
    s = (len(exposure_index),20)
    co2_interp =  zeros( s, dtype='f' )
    co2_ak_corrected =  zeros( s, dtype='f' )
    xCO2_true = zeros( len(exposure_index), dtype='f' )
    xCO2_true_AK = zeros( len(exposure_index), dtype='f' )
    # loop over all entries
    for i in range(len(exposure_index)):
       # print  p[i,:n_fine[i]],gases[i,3,:n_fine[i]]/gases[i,1,:n_fine[i]]
        co2_interp[i,:n_coarse[i]] = interp(p_coarse[i,:n_coarse[i]],p[i,:n_fine[i]],old_div(gases[i,3,:n_fine[i]],gases[i,1,:n_fine[i]]))
        AK_co2 = AK_all[i,:n_coarse[i],:n_coarse[i]]
        co2_ak_corrected[i,:n_coarse[i]] = co2_prior[i,:n_coarse[i]]+inner(AK_co2,co2_interp[i,:n_coarse[i]]-co2_prior[i,:n_coarse[i]])
        xCO2_true[i] = dot(h[i,:n_coarse[i]].transpose(),co2_interp[i,:n_coarse[i]])*1e6
        xCO2_true_AK[i] = dot(h[i,:n_coarse[i]].transpose(),co2_ak_corrected[i,:n_coarse[i]])*1e6

    l2_obj["truth/co2_profile/co2_truth"] = co2_interp
    l2_obj["truth/co2_profile/xco2_truth_calc"] = xCO2_true
    l2_obj["truth/co2_profile/xco2_truth_calc_AK_applied"] = xCO2_true_AK
    l2_obj["truth/co2_profile/co2_truth_AK_applied"] = co2_ak_corrected
    
    ###################################################################################
    # Use the aerosol optical depth profiles from the opt files.
    opt = h5py.File(options.opt_file, "r")
    aod = opt["Simulation/Aerosol/optical_thickness"][:][exposure_index,:,:,:]
    species = opt["Simulation/Aerosol/species_id"][:][exposure_index,:]
    num_aero = opt["Simulation/Aerosol/num_species"][:][exposure_index]
    p_hr =  l2_obj["truth/simulation/Thermodynamic/pressure_level"][:]
    n_hr = sim["Simulation/Thermodynamic/num_layers"][:][exposure_index]
    a = aod.shape
    
    od_total =  zeros((a[0],a[2],a[3]) , dtype='f' )
    od_ice =  zeros((a[0],a[2],a[3]) , dtype='f' )
    od_water =  zeros((a[0],a[2],a[3]) , dtype='f' )
    od_aerosol =  zeros((a[0],a[2],a[3]) , dtype='f' )

    od_total_lr =  zeros((a[0],a[2],19) , dtype='f' )
    od_ice_lr =  zeros((a[0],a[2],19) , dtype='f' )
    od_water_lr =  zeros((a[0],a[2],19) , dtype='f' )
    od_aerosol_lr =  zeros((a[0],a[2],19) , dtype='f' )

    p_hr2 = old_div((p_hr[:,1:]+p_hr[:,0:-1]),2)
    
    for i in range(a[0]):
        od_total[i,:,:] = aod[i,:,:,:].sum(axis=0)
        
        ind_water = zeros(0 , dtype='i' )
        ind_ice = zeros(0 , dtype='i' )
        ind_aero = zeros(0 , dtype='i' )
         # Look for different aerosol types and find indices...
        for typ_id in range(num_aero[i]):
         #   print species[i,typ_id][0:3]
            if species[i,typ_id][0:3]=='wat':
                ind_water = append(ind_water,typ_id)
            elif species[i,typ_id][0:3]=='ice':
                ind_ice = append(ind_ice,typ_id)
            else:
                ind_aero = append(ind_aero,typ_id)
        # Look for different profiles to be added into retrieval layers...
        for press_id in range(n_level[i]-1):
            ind = ((p_hr2[i,:]>=press_level[i,press_id])&(p_hr2[i,:]<press_level[i,press_id+1])).nonzero()
           
            for spec in range(a[2]):
                od_total_lr[i,spec,press_id] = aod[i,:,spec,ind].sum()
                if len(ind_water)>0:
                    for ind_a in ind_water:
                        od_water_lr[i,spec,press_id] += aod[i,ind_a,spec,ind].sum()
                if len(ind_aero)>0:
                    for ind_a in ind_aero:
                        od_aerosol_lr[i,spec,press_id] += aod[i,ind_a,spec,ind].sum()
                if len(ind_ice)>0:
                    for ind_a in ind_ice:
                        od_ice_lr[i,spec,press_id] = aod[i,ind_a,spec,ind].sum()
                
        od_water[i,:,:] = aod[i,ind_water,:,:].sum(axis=0)
        od_ice[i,:,:] = aod[i,ind_ice,:,:].sum(axis=0)
        od_aerosol[i,:,:] = aod[i,ind_aero,:,:].sum(axis=0)
     #   print "1: ", od_aerosol[i,0,:].sum(), od_aerosol[i,1,:].sum(), l2_obj["truth/log/tau_aerosol_1"][i]
     #   print "3: ", od_aerosol[i,4,:].sum(), od_aerosol_lr[i,4,:].sum(), l2_obj["truth/log/tau_aerosol_3"][i]
       
        
    l2_obj["truth/simulation/Aerosol/optical_thickness"] = aod
    l2_obj["truth/simulation/Aerosol/optical_thickness_all"] = od_total
    l2_obj["truth/simulation/Aerosol/optical_thickness_water"] = od_water
    l2_obj["truth/simulation/Aerosol/optical_thickness_ice"] = od_ice
    l2_obj["truth/simulation/Aerosol/optical_thickness_aero"] = od_aerosol
  
    l2_obj["truth/simulation/Aerosol/optical_thickness_all_lowRes"] = od_total_lr
    l2_obj["truth/simulation/Aerosol/optical_thickness_water_lowRes"] = od_water_lr
    l2_obj["truth/simulation/Aerosol/optical_thickness_ice_lowRes"] = od_ice_lr
    l2_obj["truth/simulation/Aerosol/optical_thickness_aero_lowRes"] = od_aerosol_lr
    
    ###################################################################################
    # Do the cloud screen stuff.
    if options.cloud_file == None:
        print("will skip the A-band cloud screen results")
    else:
        l2_obj.create_group("ABandCloudScreen")
        cld = h5py.File(options.cloud_file, "r")
        aband = cld["ABandCloudScreen"]
        for data in list(aband.items()):
            if len(data[1].shape)==1:
                l2_obj['ABandCloudScreen/'+data[0]]=data[1][:][exposure_index]
            elif len(data[1].shape)==2:
                l2_obj['ABandCloudScreen/'+data[0]]=data[1][:][exposure_index,:]
        
        cld.close()
    opt.close()
    sim.close()
    l2_obj.close()
    l1b_obj.close()
    print('done...')

def standalone_main():
    parser = OptionParser(usage="usage: %prog --l2 l2_spliced.h5 --l1b l1b.h5 --logFile l1b.log --sceneFile scene_xxx.log --detailedLog scene_xxx.hdf ")

    parser.add_option( "--l2", dest="l2_file",
                       metavar="FILE",
                       help="spliced L2 dataset (or single sounding, whatever you like")
    parser.add_option( "--l1b", dest="l1b_file",
                       metavar="FILE",
                       help="Name of the l1b file used in this run")
    parser.add_option( "--logFile", dest="log_file",
                       metavar="FILE",
                       help="Name of the log file for the l1b file used in this run")
    parser.add_option( "--sceneFile", dest="scene_file",
                       metavar="FILE",
                       help="Name of the scene file for the l1b file used in this run")
    parser.add_option( "--detailedLog", dest="logHDF_file",
                       metavar="FILE",
                       help="Name of the detailed log file (HDF5 format) for the l1b file used in this run")
    parser.add_option( "--optFile", dest="opt_file",
                       metavar="FILE",
                       help="Name of the *.opt file (HDF5 format) containing the aerosol optical depths profiles")
    parser.add_option( "--cloudFile", dest="cloud_file",
                       metavar="FILE",
                       help="Name of the cloud l2 file (HDF5 format) containing the cloud screeening fit")

    # Parse command line arguments
    (options, args) = parser.parse_args()
    
    if options.l2_file == None:
        parser.error('L2 file has to be specifiec')
    if options.l1b_file == None:
        parser.error('L1b file has to be specified')
    if options.log_file == None:
        parser.error('log file has to be specified')
    if options.scene_file == None:
        parser.error('scene file has to be specified')
    if options.logHDF_file == None:
        parser.error('log HDF file has to be specified')
    if options.opt_file == None:
        parser.error('opt HDF file has to be specified')
    if options.cloud_file == None:
        print("Cloud file not provided, ignoring for the moment...")
        
        
    mergeFiles(options)
        
    
if __name__ == "__main__":
    standalone_main()

#no guaranty, Christian Frankenberg
