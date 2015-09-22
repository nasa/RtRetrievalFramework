#!/usr/bin/env python
# This utility takes a file containing a list of L1B/ECMWF/Cloud filenames as well as a file
# containing a list of sounding ids. The program will then extract the data for those sounding
# ids from the various files and splice them into a new L1B and and additional files.

import os
import re
import sys
import copy
from optparse import OptionParser
from operator import itemgetter
from collections import Counter, deque
from contextlib import closing
import logging
import traceback
import subprocess
import logging
import time
import numpy
import h5py

import full_physics.acos_file as acos_file
from full_physics import write_progress
from full_physics.fill_value import FILL_VALUE

from full_physics.splice_tool_mapping import aggregator_dataset_mapping, aggregator_dataset_dest_names

DEFAULT_FN_TMPL = "{input_name}_spliced.h5"

logger = logging.getLogger()

# These ID dimensions are used by IdDatasetFinder for use by the splice tool and l2_analysis
ID_DIMENSIONS = { 'SoundingGeometry/sounding_id':  ('Frame', 'Sounding'),
                  'SoundingHeader/sounding_id':    ('Exposure',), }

class InputContainerChooser(object):
    """Returns the appropriate h5py wrapper class that can return id shape names"""

    def __new__(cls, filename, mode=None, single_id_dim=False):
        try:  
            obj = acos_file.IdDatasetFinder(filename, mode, single_id_dim=single_id_dim, srch_dim_ids=ID_DIMENSIONS)
        except acos_file.MissingDataset:
            obj = acos_file.SoundingFirstFile(filename, mode, single_id_dim=single_id_dim)

        return obj

class Timer:    
    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start
   
def find_sounding_indexes(in_filenames, in_srch_sounding_ids=None, multi_source_mode=False):
    """Returns a list of ids for each passed object as well as a list of indexes for each file
    where these sounding ids come from the list passed as an argument"""
    
    ids_by_file = []
    indexes_by_file = []
    all_sounding_ids = set()
    
    if in_srch_sounding_ids != None:
        in_srch_sounding_ids = set(numpy.array(in_srch_sounding_ids, dtype=numpy.int64))

    for file_idx, curr_file in enumerate(in_filenames):
        # Write progress at beginning of loop so code doesn't look stalled out
        # -1 since we have not completed work on this iteration
        write_progress(file_idx-1, len(in_filenames))
            
        with closing(InputContainerChooser(curr_file, single_id_dim=not multi_source_mode)) as hdf_obj:
            file_sounding_ids = set(numpy.ravel(hdf_obj.get_sounding_ids()))

            if in_srch_sounding_ids == None:
                matching_ids = file_sounding_ids
            else:
                matching_ids = file_sounding_ids.intersection(in_srch_sounding_ids)
            all_sounding_ids.update(matching_ids)

            matching_ids = sorted(matching_ids)
            matching_indexes = [ hdf_obj.get_sounding_indexes(curr_id) for curr_id in matching_ids ] 

            ids_by_file.append( tuple(matching_ids) )
            indexes_by_file.append( tuple(matching_indexes) )

        write_progress(file_idx, len(in_filenames))

    total_found = len(all_sounding_ids)
    return (ids_by_file, indexes_by_file, total_found)

class DatasetInformation(object):

    def __init__(self, dataset_name, rename_mapping=False):
        # Strip leading and trailing slashes just in case
        self.inp_name = dataset_name.strip("/")

        base_name = self.inp_name.split("/")[-1]
        if rename_mapping and aggregator_dataset_mapping.has_key(base_name):
            self.out_name = aggregator_dataset_mapping[base_name]["name"]
        else:
            self.out_name = self.inp_name

        # Default values
        self.id_names = None
        self.inp_shape_names = None
        self.out_shape_names = None
        self.inp_id_dims = None

        # Store shapes for each file so we can attach this as an attribute
        self.inp_file_shapes = []

        # Which h5py objects contain the dataset
        self.inp_filenames = []

        # Type and shape for output dataset
        # out_shape should be the max of all file shapes in all input file
        self.out_type = None
        self.out_shape = None 

    def __str__(self):
        res = ""

        res += "inp_name: %s\n" % self.inp_name
        res += "out_name: %s\n" % self.out_name
        res += "id_names: %s\n" % self.id_names
        res += "inp_shape_names: %s\n" % self.inp_shape_names
        res += "out_shape_names: %s\n" % self.out_shape_names
        res += "inp_id_dims: %s\n" % self.inp_id_dims
        res += "inp_file_shapes: %s\n" % self.inp_file_shapes
        res += "inp_filenames: %d files\n" % len(self.inp_filenames)
        res += "out_type: %s\n" % self.out_type
        res += "out_shape: %s\n" % self.out_shape

        return res

    def add_instance(self, file_obj, dataset_obj):
        "Add an instance of dataset from new file"

        # Track the file objects considered for this dataset
        self.inp_filenames.append(file_obj.filename)

        if self.inp_shape_names == None: 
            # Find the names of dimensions
            self._find_dim_names(file_obj, dataset_obj)

            # Find where id dimensions are indexed
            self._find_id_dims(file_obj, dataset_obj)

        self._determine_output_shape(file_obj, dataset_obj)

    def _find_dim_names(self, file_obj, dataset_obj):
        "Finds the dimension names of the dataset that should be used in the output"

        if self.inp_name.find("Iteration/") >= 0:
            # Set up Iteration/ datasets based on their non Iteration/ versions's shape
            # Since Iteration/ datasets do not have a shape of their own
            non_iter_name = self.inp_name.replace("Iteration/", "")
            self.inp_shape_names = file_obj.get_data_shape(non_iter_name)
        else:
            self.inp_shape_names = file_obj.get_data_shape(dataset_obj)

    def _find_id_dims(self, file_obj, dataset_obj):
        "Finds where in the dataset's dimensions the id dimensions are indexed"

        # Name of dimensions denoting id of datasets, eg sounding_id
        self.id_names = file_obj.get_id_dim_names()

        # Retrieve directly from file as self.inp_shape_names is
        # the output file names, might not match source files
        ds_dim_names = file_obj.get_data_shape(dataset_obj)

        self.inp_id_dims = []
        for curr_id_name in self.id_names:
            if curr_id_name in ds_dim_names:
                self.inp_id_dims.append(ds_dim_names.index(curr_id_name))

    def _determine_output_shape(self, file_obj, dataset_obj):
        # Try and locate files that contain the dataset
        self.inp_file_shapes.append( dataset_obj.shape ) 

        # Only consider shape of additional items if the dataset has id dimensions and is not just a file we will
        # copy from the first available
        if self.out_shape != None and len(self.inp_id_dims) > 0:
            if dataset_obj.dtype.kind == 'S' and self.out_type.kind == 'S':
                # Use the largest string type found
                if dataset_obj.dtype.itemsize > self.out_type.itemsize:
                    self.out_type = dataset_obj.dtype
            elif dataset_obj.dtype.kind != self.out_type.kind:
               raise TypeError("%s does not have the same datatype as all other files for dataset: %s." % (file_obj.filename, self.inp_name))

            #  Find maximum of all datasets so far
            self.out_shape = [ max(self.out_shape[i], curr_shape) for i, curr_shape in enumerate(dataset_obj.shape) ]

        elif self.out_shape != None and len(self.inp_id_dims) == 0:
            # Already have a shape and this is not an id dimension so just ignore any shape updates
            return
        else:
            # Assign the type of the data an make sure all input file use the same one
            self.out_type = dataset_obj.dtype
            self.out_shape = list(dataset_obj.shape)

    def input_data_indexes(self, dataset_obj, indexes):
        "Faster version of what is in acos_file.SoundingDataFile"

        if self.inp_shape_names == None:
            raise ValueError('No shape names are defined for dataset: %s in file: %s' % (dataset_name, self.filename))
        elif len(self.inp_shape_names) != len(dataset_obj.shape):
            raise ValueError('Length of shape names: %s does not match length of data: %s for: %s' % (len(self.inp_shape_names), len(dataset_obj.shape), dataset_obj.name))

        if len(self.id_names) > 1:
            index_get = [ itemgetter(idx) for idx in range(len(self.id_names)) ]
        else:
            index_get = [ lambda f:f ]

        data_dim_indexes = []
        for dim_idx, shp_name in enumerate(self.inp_shape_names):
            if shp_name in self.id_names:
                shape_idxs = tuple(map(index_get[self.id_names.index(shp_name)], indexes))

                if len(shape_idxs) == 1:
                    shape_idxs = shape_idxs[0]
                data_dim_indexes.append(shape_idxs)
            else:
                data_dim_indexes.append( slice(dataset_obj.shape[dim_idx]) )

        # Eliminate duplicate shape names
        for name, num_occurance in Counter(self.inp_shape_names).items():
            if num_occurance > 1:
                for name_count in range(1, num_occurance+1):
                    self.inp_shape_names[self.inp_shape_names.index(name)] = "%s_%d" % (name, name_count)

        return tuple(data_dim_indexes)
 
    def output_dataset_shape(self, splice_size=None, collapse_id_dim=False):
        # Make a copy of dataset shape size and names as they will be modified
        dst_shape = self.out_shape 
        self.out_shape_names = list(self.inp_shape_names)
        splice_dim = self.inp_id_dims

        if splice_dim != None and len(splice_dim) > 0 and splice_size != None:
            # Set first splice dim to splice size
            dst_shape[splice_dim[0]] = splice_size

            # Concatenate input shape names for id dims into one name
            # Go through and collapse splice_dim dim(s) into a single one of
            # sized according to splice_dim
            if collapse_id_dim:
                if self.out_shape_names != None:
                    out_id_name = ''.join([ self.out_shape_names[didx] for didx in splice_dim ])
                    self.out_shape_names[splice_dim[0]] = out_id_name

                # Remove any other splice dims,
                # Iterate backwards since we are
                # popping from dst_shape
                for s_dim in splice_dim[1:][::-1]:
                    dst_shape.pop(s_dim)
                    if self.out_shape_names != None:
                        self.out_shape_names.pop(s_dim)

        # Make sure any zero sized datasets are indicated as unlimited
        # by identifying them as None in the maxshape argument
        max_shape = map(lambda s: None if s == 0 else s, dst_shape)

        return dst_shape, max_shape

def determine_datasets_info(in_filenames, desired_datasets=None, search_all_files=True, multi_source_mode=False, rename_mapping=False):
    # Arrays returned
    copy_datasets = {}

    logger.info("Determining datasets information")

    # Always check first file 
    check_files = [ in_filenames[0] ]

    # Optionally search for all other files for datasets and sizes
    if search_all_files and len(in_filenames) > 1:
        logger.info("Will search all files for dataset information")
        check_files += in_filenames[1:]

    ignored_datasets = []
    file_id_names = []
    with Timer() as t:
        for file_idx, curr_file in enumerate(check_files):
            with closing(InputContainerChooser(curr_file, single_id_dim=not multi_source_mode)) as hdf_obj:
                def check_dataset(ds_name, ds_obj):
                    if isinstance(ds_obj, h5py.Group) or ds_name in ignored_datasets:
                        return None

                    ds_info = copy_datasets.get(ds_name, DatasetInformation(ds_name, rename_mapping))

                    # Search list of datasets based upon the destination dataset name, which is the same ass
                    # the input name when rename mapping is not used
                    if desired_datasets != None and len(filter(lambda ds: re.search(ds_info.out_name, ds), desired_datasets)) == 0:
                        logger.debug('Ignoring dataset not in list of datasets to consider: "%s"' % ds_name)
                        ignored_datasets.append(ds_name)
                        return None

                    ds_info.add_instance(hdf_obj, ds_obj)
                    if len(ds_info.inp_shape_names) > 0:
                        copy_datasets[ds_name] = ds_info
                    else:
                        logger.warning( "Ignoring dataset which is not an array: %s" % ds_name )

                hdf_obj.visititems(check_dataset)

                if len(check_files) > 1:
                    write_progress(file_idx, len(check_files))

                # Add file id names, for checking for the type of the file, w/o opening it again
                file_id_names.append(hdf_obj.get_id_dim_names())

    logger.info("Analyzed %d files and found %d unique datasets in %.03f seconds" % (len(check_files), len(copy_datasets), t.interval))

    return (copy_datasets, file_id_names)

def create_output_dataset(out_hdf_obj, dataset_info, splice_size=None, collapse_id_dim=False):
    """Duplicates a dataset from the input file into the output hdf object as it exists
    except for its dimensions"""

    logger.debug("Creating new output dataset: %s" % dataset_info.out_name)

    dst_shape, max_shape = dataset_info.output_dataset_shape(splice_size, collapse_id_dim) 

    # Split name into two and create a group if needed
    # then create the desire dataset or load existing group
    ds_name_clean = dataset_info.out_name.lstrip('/').rstrip('/')
    if ds_name_clean.find("/") > 0:
        dst_group, dst_name = dataset_info.out_name.lstrip('/').rstrip('/').split('/', 1)
        out_group_obj = out_hdf_obj.require_group(dst_group) 
    else:
        dst_group = ""
        dst_name = ds_name_clean
        out_group_obj = out_hdf_obj

    # Fill new dataset with the correct fill value based on type
    if dataset_info.out_type != numpy.object and dataset_info.out_type.type != numpy.string_:
        fill_type = dataset_info.out_type.type

        if FILL_VALUE.has_key(fill_type):
            dataset_fill = FILL_VALUE[fill_type]
        else:
            logger.warning("Could not find specific fill value for dataset: %s of type %s" % (dst_name, fill_type))
            dataset_fill = None
    else:
        # Use default fill for string types
        fill_type = None
        dataset_fill = None

    logger.debug( "Creating new dataset: %s/%s sized: %s with fill type: %s and value: %s" % (dst_group, dst_name, dst_shape, fill_type, dataset_fill) )
    try:
        out_dataset_obj = out_group_obj.create_dataset(dst_name, data=numpy.empty(dst_shape, dtype=dataset_info.out_type), maxshape=max_shape, compression="gzip", compression_opts=2, fillvalue=dataset_fill)
    except RuntimeError as exc:
        raise RuntimeError("Error creating dataset %s/%s: %s" % (dst_group/dst_name, exc))

    # Now create copied attributes from original dataset
    # Just copy from first for now, leave code to do multiple if needed
    for curr_file in dataset_info.inp_filenames[0:1]:
        with closing(h5py.File(curr_file, 'r')) as curr_hdf_obj:
            curr_dataset_obj = curr_hdf_obj[dataset_info.inp_name]
            for attr_name, attr_value in curr_dataset_obj.attrs.items():
                # Skip if copied already from a file, assuming all files
                # have same attributes for now, we just will get
                # all uniquely named ones
                if attr_name in out_dataset_obj.attrs.keys():
                    continue
        
                # If the dtype of the attribute is an object dtype, then assume
                # its a variable length string
                if hasattr(attr_value, "dtype") and attr_value.dtype.kind == "O":
                    logger.debug('Copying variable length string attribute: "%s" with value: "%s"' % (attr_name, attr_value[0]))
                    vlen_dt = h5py.special_dtype(vlen=str)
                    out_dataset_obj.attrs.create(attr_name, attr_value, dtype=vlen_dt)
                else:
                    logger.debug('Copying attribute: "%s" with value: "%s"' % (attr_name, attr_value))
                    out_dataset_obj.attrs.create(attr_name, attr_value)

    # Add extra information for dataset, overwrite an existing shape, because we may have
    # reshaped the data
    if dataset_info.out_shape_names:
        out_dataset_obj.attrs["Shape"] = numpy.array(["_".join(dataset_info.out_shape_names) + "_Array"]) 

    # Only add dimensions for individual files if files had differing dimension sizes
    file_shapes = dataset_info.inp_file_shapes
    if len(file_shapes) > 0 and (len(filter(lambda x: x == file_shapes[0], file_shapes)) != len(file_shapes)):
        try:
            file_shapes_arr = numpy.array(file_shapes) 
            out_dataset_obj.attrs.create("File_Dimensions", data=file_shapes_arr)
        except RuntimeError:
            # This will happen if the dataset would be too big
            logger.warning("Unable to add File_Dimensions attribute to %s dataset" % dataset_info.out_name)

    return out_dataset_obj

def get_dataset_slice(acos_hdf_obj, in_dataset_obj, dataset_info, in_data_idx, out_shape):
    """Copys dataset values from one dataset object to another, but only certain indexes along a
    specific dimension of the data"""

    # Determine how to extact data other than the splice dimension
    in_dataset_indexes = dataset_info.input_data_indexes(in_dataset_obj, in_data_idx)

    # Obtain selected data for copying into output dataset
    try:
        if len(in_dataset_indexes) == 1 and not isinstance(in_dataset_indexes[0], slice):
            in_data = in_dataset_obj[:][numpy.array(in_dataset_indexes[0])]
        else:
            in_data = in_dataset_obj[:][tuple(in_dataset_indexes)]
    except IOError as exc:
        raise IOError("Can not read dataset %s from file %s: %s" % (dataset_info.inp_name, acos_hdf_obj.filename, exc))

    # Set sliced data into output dataset
    if numpy.product(in_data.shape) > numpy.product(out_shape):
        logger.warning("Dataset %s requires destructive resizing" % (dataset_info.out_name))
        logger.debug("At indexes %s resizing source data of shape %s to %s." % (in_data_idx, in_data.shape, out_shape)) 
        stored_data = numpy.resize(in_data, out_shape)
    else:
        stored_data = in_data.reshape(out_shape)

    return stored_data

def create_dest_file_datasets(out_hdf_obj, in_filenames, inp_datasets_info, out_num_soundings, copy_all=False, multi_source_mode=False):
    logger.info("Creating output file datasets")

    # Loop over dataset names/shapes
    # Create datasets and copy any non id_shape datasets
    num_datasets = len(inp_datasets_info.keys())
    out_dataset_objs = {}
    for curr_index, (curr_dataset_name, curr_dataset_info) in enumerate(inp_datasets_info.items()):
        write_progress(curr_index-1, num_datasets)

        id_dimension = curr_dataset_info.inp_id_dims

        if out_hdf_obj.get(curr_dataset_info.out_name, None) != None:
            logger.debug("Destination dataset %s for %s already exists in output file, ignoring duplicate" % (curr_dataset_info.out_name, curr_dataset_name))
            del inp_datasets_info[curr_dataset_name]

        elif len(id_dimension) > 0: 
            # Create a dataset that does have an id dimension
            out_dataset_objs[curr_dataset_name] = create_output_dataset(out_hdf_obj, curr_dataset_info, out_num_soundings, collapse_id_dim=multi_source_mode)
        
        elif copy_all:
            logger.debug("Non sounding id dataset, copying from first file: %s" % curr_dataset_name)
                
            # Copy dataset with no sounding dimension directly from first file
            out_dataset_obj = create_output_dataset(out_hdf_obj, curr_dataset_info, collapse_id_dim=multi_source_mode)
           
            # Search for first file that has the desired dataset. Can't just take
            # from first in cases when splicing dissimilar files
            for inp_fn in in_filenames:
                with closing(h5py.File(inp_fn, 'r')) as hdf_obj:
                    if hdf_obj.get(curr_dataset_name, None):
                        # Do not try and copy empty datasets
                        if hdf_obj[curr_dataset_name].len() > 0:
                            out_dataset_obj[:] = hdf_obj[curr_dataset_name][:]
                        inp_datasets_info.pop(curr_dataset_name)
                        break
        else:
            del inp_datasets_info[curr_dataset_name]
            logger.debug( "Ignoring non sounding id dataset: %s" % curr_dataset_name )

        write_progress(curr_index, num_datasets)
    
    return out_dataset_objs

def output_indexes_shape(in_dataset_obj, dataset_info, out_data_idx, collapse_id_dim=False):
    # Create output dataset indexes, removing any concatenated dimensions 
    splice_dim = dataset_info.inp_id_dims 
    out_dataset_idx = [ slice(dim_size) for dim_size in in_dataset_obj.shape ]
    out_dataset_idx[splice_dim[0]] = out_data_idx
    if collapse_id_dim:
        for sidx in splice_dim[1:][::-1]:
            out_dataset_idx.pop(sidx)

    # Get destination dataset shape so we can reshape the incoming data to match
    out_shape = []
    for idx in out_dataset_idx:
        if type(idx) is slice:
            start = idx.start and idx.start or 0
            stop = idx.stop
            if idx.step:
                raise NotImplementedError("stepped slices not supported yet")
            out_shape.append(stop - start)
        elif type(idx) is int:
            out_shape.append(1)
        else:
            raise ValueError("Unknown type %s in indexes: %s" % (type(idx), out_dataset_idx))

    return out_dataset_idx, out_shape

def get_datasets_for_file(curr_file, curr_snd_indexes, inp_datasets_info, output_index, multi_source_mode=False):
    # Ignore current file if no soundings
    if len(curr_snd_indexes) == 0:
        return     

    with closing(InputContainerChooser(curr_file, single_id_dim=not multi_source_mode)) as curr_hdf_obj:
        # Copy each dataset for current file
        for curr_dataset_name, curr_dataset_info in inp_datasets_info.items():

            curr_ds_obj = curr_hdf_obj.get(curr_dataset_name, None)
            if curr_ds_obj != None:
                output_slice = slice(output_index, output_index+len(curr_snd_indexes))

                # Only actually perform the copy if the object does not contain an empty dimension
                # otherwise we leave the created output dataset empty
                if not 0 in curr_ds_obj.shape:
                    # Copy dataset's relevant contents into output
                    try:
                        out_dataset_idx, out_shape = output_indexes_shape(curr_ds_obj, curr_dataset_info, output_slice, collapse_id_dim=multi_source_mode)

                        stored_data = get_dataset_slice(curr_hdf_obj, curr_ds_obj, curr_dataset_info, curr_snd_indexes, out_shape)
                        yield curr_dataset_name, out_dataset_idx, stored_data
      
                    except (ValueError, IOError) as e:
                        logger.error('Error copying dataset: "%s", dataset may be corrupt: %s' % (curr_dataset_name, e))
                        logger.debug('Exception: %s' % traceback.format_exc())
                else:
                    logger.debug('Dataset: "%s" has an empty dimension in its shape: %s, not actually copying any values from: %s' % (curr_dataset_name, curr_ds_obj.shape, curr_file))

            else: 
                # Ok if not found, just treat as empty dataset
                logger.debug('Dataset: "%s" is missing, not actually copy any values from: %s' % (curr_dataset_name, curr_file))

def copy_multiple_datasets(out_dataset_objs, in_filenames, inp_snd_indexes, file_id_names, inp_datasets_info, multi_source_mode=False):
    """Given an input and output hdf object, copies the specified soundings ids for the specified
    dataset names along the specified shape dimention to the output file"""

    logger.info("Copying data to output")
 
    # Loop over dataset names/shapes
    output_index = 0
    for file_index, (curr_file, curr_snd_indexes) in enumerate(zip(in_filenames, inp_snd_indexes)):
        logger.debug( 'Copying %d sounding(s) from %s (%d / %d)' % (len(curr_snd_indexes), os.path.basename(curr_file), file_index+1, len(in_filenames)) )

        write_progress(file_index-1, len(in_filenames))

        for curr_dataset_name, out_dataset_idx, stored_data in get_datasets_for_file(curr_file, curr_snd_indexes, inp_datasets_info, output_index, multi_source_mode):
            out_dataset_objs[curr_dataset_name].__setitem__(tuple(out_dataset_idx), stored_data)

        write_progress(file_index, len(in_filenames))

        # Increment where to start in the index into output sounding indexes
        # Only increment if not in multi-source mode because when using different
        # inputs to a single output we will be pulling the same sounding ids (hopefully)
        if not multi_source_mode:
            output_index += len(curr_snd_indexes) 

def determine_multi_source(check_files):

    # Looks through input hdf files and see if we are in multi-source combination mode, ie different types of files combined into one
    multi_source_mode = False

    file_id_names = []
    for curr_file in check_files:
        with closing(InputContainerChooser(curr_file)) as hdf_obj:
            file_id_names.append( hdf_obj.get_id_dim_names() )

    # Lets make sure the inputs are either all different types
    # or all the same, otherwise if some are the same,
    # bad things could happen
    id_match_counts = numpy.array([ len(filter(lambda x: x == idn, file_id_names)) for idn in file_id_names ])

    if(all(id_match_counts == len(id_match_counts))):
        logger.debug("Not using multi source mode")
    elif(len(id_match_counts) > 1 and all(id_match_counts >= 1) and all(id_match_counts < len(id_match_counts))):
        multi_source_mode = True
        logger.debug("Using multi source mode")
    else:
        raise Exception("Error determining multi source mode with id match counts: %s" % id_match_counts)

    return multi_source_mode

def splice_acos_hdf_files(id_src_files, sounding_ids, output_id_src, output_addl={}, splice_all=False, desired_datasets_list=None, search_all_files=True, rename_mapping=False):

    # Determine if input files are of different types
    multi_source_mode = determine_multi_source(id_src_files)
    
    for addl_files in output_addl.values():
        if len(addl_files) != len(id_src_files):
            raise Exception("Additional files lists must have same number of files as l1b file list")

    dest_id_obj = h5py.File(output_id_src, 'w')

    copy_filenames_dict = {}
    copy_filenames_dict[dest_id_obj] = id_src_files
    
    for addl_outname, addl_file_list in output_addl.items():
        try:
            dest_addl_obj = h5py.File(addl_outname, 'w')
        except IOError as e:
            raise IOError("Could not create output file: %s\n%s" % (addl_outname, e))
        copy_filenames_dict[dest_addl_obj] = addl_file_list

    # Match sounding ids to files and indexes using L1B files
    logger.info("Finding sounding indexes from source file sounding ids")
    with Timer() as t:
        (sounding_id_matches, sounding_index_matches, total_num_soundings) = find_sounding_indexes(id_src_files, sounding_ids, multi_source_mode=multi_source_mode)
    logger.info("Found %d total soundings from %d files in %.03f seconds" % (total_num_soundings, len(id_src_files), t.interval))

    # Sort filename based lists
    # Pack together, sort, then unpack
    for dest_obj in copy_filenames_dict.keys():
        sort_lists = zip(sounding_id_matches, sounding_index_matches, copy_filenames_dict[dest_obj])
        sort_lists.sort(key=itemgetter(0)) # Sort by sounding id
        copy_filenames_dict[dest_obj] = map(itemgetter(2), sort_lists)

    # Do not extract these till after we have sorted all dest_objs
    # In the loop above these two will get sorted multiple times,
    # Serving as the key for the sort of the other objs. But we
    # have to leave the input for sort unsorted to get the same
    # sorting behavior for each dest obj
    sounding_id_matches = map(itemgetter(0), sort_lists)
    sounding_index_matches = map(itemgetter(1), sort_lists)
       
    # Copy over all datasets in files
    for dest_obj, in_filenames in copy_filenames_dict.items():
        logger.info( "Splicing into: %s" % dest_obj.filename )

        # Trim out files that have no matched soundings
        if len(in_filenames) > 1:
            where_has_snd_ids = filter(lambda i: len(sounding_index_matches[i]) > 0, range(len(sounding_index_matches)))
            if len(where_has_snd_ids) == 0:
                raise Exception("Failed to find any files with matching sounding ids in them")
            get_valid = itemgetter(*where_has_snd_ids)
            in_filenames = get_valid(in_filenames)
            copy_indexes = get_valid(sounding_index_matches)

            # Make a list if getter found only one item and made it just a single string
            if not hasattr(in_filenames, "__iter__"):
                in_filenames = [ in_filenames ]
                copy_indexes = [ copy_indexes ]
        else:
            copy_indexes = sounding_index_matches
        
        # Determine information on datasets found in input files
        (copy_datasets, file_id_names) = determine_datasets_info(in_filenames, desired_datasets_list, search_all_files, multi_source_mode=multi_source_mode, rename_mapping=rename_mapping)

        # Create datasets in the output file so we have somewhere to dump the input
        out_dataset_objs = create_dest_file_datasets(dest_obj, in_filenames, copy_datasets, total_num_soundings, splice_all, multi_source_mode)
        dest_obj.flush()

        # Now copy desired datasets from source file to destination files
        with Timer() as t:
            copy_multiple_datasets(out_dataset_objs, in_filenames, copy_indexes, file_id_names, copy_datasets, multi_source_mode=multi_source_mode)

        logger.info("Datasets copying took %.03f seconds" % (t.interval))
        dest_obj.close()

def common_default_filename(fullnames):
    basenames = [os.path.basename(fn) for fn in fullnames]
    common_name = os.path.commonprefix(basenames)
    common_name = common_name.rstrip("_").rstrip("-") # convience to remove maybe a common _ mark
    return DEFAULT_FN_TMPL.format(input_name=os.path.splitext(common_name)[0])
    
def standalone_main():
    parser = OptionParser(usage="usage: %prog [ -i <id_file_plus_addl_list_file> | <filename_1> <filename_2> ... ] [ -s <sounding_id_file> ] [ -o <first_output_file> ] [ -a <additional_output_name> ... ]")

    parser.add_option( "-i", "--input_files_list", dest="input_files_list",
                       metavar="FILE",
                       help="text file with main input filename plus any additional filenames to splice separated by whitespace on each line")

    parser.add_option( "-s", "--sounding_id_file", dest="sounding_id_file",
                       metavar="FILE",
                       help="file containing list of soundings for destination file(s)")

    parser.add_option( "-o", "--output_file", dest="first_output_file",
                       metavar="FILE",
                       help="name for main output filename, aka the files containing the sounding ids")

    parser.add_option( "-a", "--addl_output_file", dest="addl_output_file",
                       metavar="FILE",
                       action="append",
                       help="spliced output filenames other than main output file, specify one for each additional file in input file list")

    parser.add_option( "-d", "--datasets_list_file", dest="datasets_list_file",
                       metavar="FILE",
                       help="file containing list of only datasets to consider for copying. If rename_mapping is enabled then the names are matched on their destination dataset name")

    parser.add_option( "-r", "--rename_mapping", dest="rename_mapping",
                       action="store_true",
                       default=False,
                       help="rename datasets into output file according to internal mapping table as they would appear in the L2Agg PGE")

    parser.add_option( "--agg_names_filter", dest="agg_names_filter",
                       action="store_true",
                       default=False,
                       help="include only dataset names that would appear in the L2Agg PGE. Its only makes sense to use this option with --rename_mapping")

    parser.add_option( "--no-search-all", dest="search_all_files",
                       action="store_false",
                       default=True,
                       help="search all files to discover all possible datasets and maximum shapes")

    parser.add_option( "--splice-all", dest="splice_all",
                       action="store_true",
                       default=False,
                       help="splice all datasets, including those which do not have a sounding dimension. Note that datasets without an explicit handler and no sounding dimension are simply copied from the first file.")

    parser.add_option( "-v", "--verbose", dest="verbose",
                       action="store_true",
                       default=False,
                       help="enable verbose informational reporting")


    # Parse command line arguments
    (options, args) = parser.parse_args()

    if len(args) == 0 and options.input_files_list == None:
        parser.error('input list file must be specified')

    # Set up logging
    if options.verbose:
        # Include HDF5 errors in output
        h5py._errors.unsilence_errors()

        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    console = logging.StreamHandler(stream=sys.stdout)
    logger.addHandler(console)

    id_src_files = []
    addl_inp_files = []
    if len(args) > 0:
        # Make sure to only add unique filenames
        # Otherwise indexing can be messed up
        for curr_fn in args:
            real_fn = os.path.realpath(curr_fn)
            if real_fn not in id_src_files:
                id_src_files.append(real_fn)
            else:
                logger.error("Duplicate filename: %s" % real_fn)

        addl_inp_files = [ tuple() for fn in id_src_files ]
    else:
        last_seen_addl_len = 0
        line_count = 0
        with open(options.input_files_list) as input_list_fo:
            for row_string in input_list_fo:
                line_count += 1

                row_filenames = row_string.split()
                id_filename = os.path.realpath(row_filenames.pop(0))
                if not id_filename in id_src_files:
                    id_src_files.append(id_filename)
                else:
                    logger.error("Duplicate filename: %s" % id_filename)

                # All tuples in addl_inp_files have same length by virtual of this check
                if last_seen_addl_len> 0 and last_seen_addl_len != len(row_filenames):
                    parser.error("Each line of input files list must have same number of filenames. Error at line %d" % line_count)

                addl_inp_files.append( tuple(row_filenames) )
                last_seen_addl_len = len(row_filenames)

    if options.first_output_file == None:
        first_output_file = common_default_filename(id_src_files)
    else:
        first_output_file = options.first_output_file

    addl_out_names = []
    if options.addl_output_file == None:
        # Use common name among inputted files to make an output file
        for addl_idx in range(len(addl_inp_files[0])):
            fullnames = map(itemgetter(addl_idx), addl_inp_files)
            addl_out_names.append( common_default_filename(fullnames) )
                                   
    elif len(options.addl_output_file) != len(addl_inp_files[0]):
        raise parser.error("When specifying additional output file names, must specify one for each type that shows up in list file")
    else:
        addl_out_names = options.addl_output_file

    addl_out_dict = {}
    for addl_idx, addl_name in enumerate(addl_out_names):
        addl_out_dict[addl_name] = map(itemgetter(addl_idx), addl_inp_files)

    if options.sounding_id_file != None:
        sounding_ids = []
        with open(options.sounding_id_file) as id_fo:
            for snd_id_line in id_fo:
                sounding_ids.append( snd_id_line.strip() )
    else:
        logger.info("No sounding ids file supplied, aggregating all ids from all files")
        sounding_ids = None

    copy_datasets_list = None
    if options.datasets_list_file != None:
        with open(options.datasets_list_file) as ds_list_fo:
            copy_datasets_list = [ ln.strip() for ln in ds_list_fo.readlines() ]
    if options.agg_names_filter:
        if copy_datasets_list == None:
            copy_datasets_list = aggregator_dataset_dest_names
        else:
            copy_datasets_list += aggregator_dataset_dest_names

    splice_acos_hdf_files(id_src_files, sounding_ids, first_output_file, addl_out_dict, splice_all=options.splice_all, desired_datasets_list=copy_datasets_list, search_all_files=options.search_all_files, rename_mapping=options.rename_mapping)

if __name__ == "__main__":
    standalone_main()
