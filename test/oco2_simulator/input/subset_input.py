# This short program is used to subset the full orbit files to a single
# sounding. I'll record this in case it is useful again in the future,
# but the data has already been generated.
import h5py

nrow = 1
row_start = 2693
tot_row = 8344

def subset_object(obj_in, obj_out):
    row = slice(row_start,row_start + nrow)
    for f in obj_in:
        if isinstance(obj_in[f], h5py.Dataset):
            s = list(obj_in[f].shape)
            if(s[0] == tot_row):
                s[0] = nrow
                obj_out.create_dataset(f, s, obj_in[f].dtype, 
                                       compression='gzip', compression_opts = 9)
                obj_out[f][0:nrow,...] = obj_in[f][row, ...]
            else:
                obj_out.create_dataset(f, s, obj_in[f].dtype,
                                       compression='gzip', compression_opts = 9)
                obj_out[f][...] = obj_in[f][...]
        if isinstance(obj_in[f], h5py.Group):
            gout = obj_out.create_group(f)
            subset_object(obj_in[f], gout)


fin = h5py.File("OCO2_meteorology_NDa_20100909_252_r60CSUSim01a_012.hdf")
fout = h5py.File("oco2_sim_met.h5", "w")
subset_object(fin, fout)

fin = h5py.File("OCO2_sim_NDa_20100909_252_r60CSUSim01b7_012.hdf")
fout = h5py.File("oco2_sim_l1b.h5", "w")
subset_object(fin, fout)

nrow = 4
row_start = 2693 * nrow
tot_row = 8344 * nrow
fin = h5py.File("OCO2_scene_NDa_20100909_252_r60CSUSim01a_012.hdf")
fout = h5py.File("oco2_sim_scene.h5", "w")
subset_object(fin, fout)

