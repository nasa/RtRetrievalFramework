#!/usr/bin/env python
# This utility takes a file containing a list of L1B/ECMWF/Cloud filenames as well as a file
# containing a list of sounding ids. The program will then extract the data for those sounding
# ids from the various files and splice them into a new L1B and and additional files.

import os
import re
import sys
from copy import copy
from argparse import ArgumentParser
from operator import itemgetter
from collections import Counter
from contextlib import closing
import logging
import traceback
import time

import numpy
import h5py

import full_physics.acos_file as acos_file
from full_physics import write_progress
from full_physics.fill_value import FILL_VALUE

from full_physics.splice_tool_mapping import aggregator_dataset_mapping, aggregator_dataset_dest_names

DEFAULT_OUTPUT_FILENAME = "spliced_output.h5"

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

        # These don't change between files, load on the first encounter
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
        out_shape_names = list(self.inp_shape_names)
        splice_dim = self.inp_id_dims

        if splice_dim != None and len(splice_dim) > 0 and splice_size != None:
            # Set first splice dim to splice size
            dst_shape[splice_dim[0]] = splice_size

            # Concatenate input shape names for id dims into one name
            # Go through and collapse splice_dim dim(s) into a single one of
            # sized according to splice_dim
            if collapse_id_dim:
                if out_shape_names != None:
                    out_id_name = ''.join([ out_shape_names[didx] for didx in splice_dim ])
                    out_shape_names[splice_dim[0]] = out_id_name

                # Remove any other splice dims,
                # Iterate backwards since we are
                # popping from dst_shape
                for s_dim in splice_dim[1:][::-1]:
                    dst_shape.pop(s_dim)
                    if out_shape_names != None:
                        out_shape_names.pop(s_dim)

        # Make sure any zero sized datasets are indicated as unlimited
        # by identifying them as None in the maxshape argument
        max_shape = map(lambda s: None if s == 0 else s, dst_shape)

        return dst_shape, max_shape, out_shape_names

class SourceInformation(object):

    def __init__(self, src_filenames, sounding_ids=None, desired_datasets=None, rename_mapping=False, multi_source_types=False):

        # Filenames to splice from, will be modified to remove any which have no sounding ids present
        # Ensure it is a list so we can modify it as we go
        self.src_filenames = list(src_filenames)

        if sounding_ids != None:
            self.inp_sounding_ids = set(numpy.array(sounding_ids, dtype=numpy.int64))
        else:
            self.inp_sounding_ids = None

        # All sounding ids found from source files that match input list
        self.found_sounding_ids = set()

        # Flags
        self.desired_datasets = desired_datasets
        self.rename_mapping = rename_mapping
        self.multi_source_types = multi_source_types

        # Sounding ids per file and their indexes in those files
        self.file_sounding_ids = []
        self.file_indexes = []

        # Mapping of dataset names to their information object
        self.datasets_info = {}
        
        # Datasets to ignore on subsequent passes, datasets added to this list to speed up
        # processing and eliminate the need to compare against the desired dataset list on
        # each encounter
        self._ignored_datasets = []

    @property
    def num_soundings(self):
        return len(self.found_sounding_ids)

    def process_file_ids(self, hdf_obj):
        curr_sounding_ids = set(numpy.ravel(hdf_obj.get_sounding_ids()))

        if self.inp_sounding_ids == None:
            matching_ids = curr_sounding_ids
        else:
            matching_ids = sorted(curr_sounding_ids.intersection(self.inp_sounding_ids))

        if len(matching_ids) > 0:
            self.found_sounding_ids.update(matching_ids)
            
            matching_indexes = [ hdf_obj.get_sounding_indexes(curr_id) for curr_id in matching_ids ] 

            self.file_sounding_ids.append( tuple(matching_ids) )
            self.file_indexes.append( tuple(matching_indexes) )

        return len(matching_ids)

    def process_dataset_info(self, hdf_obj):
        def _check_dataset(ds_name, ds_obj):
            if isinstance(ds_obj, h5py.Group) or ds_name in self._ignored_datasets:
                return None

            ds_info = self.datasets_info.get(ds_name, None)
            if ds_info == None:
                ds_info = DatasetInformation(ds_name, self.rename_mapping)

            # Search list of datasets based upon the destination dataset name, which is the same as
            # the input name when rename mapping is not used
            if self.desired_datasets != None and len(filter(lambda ds: re.search(ds_info.out_name, ds), self.desired_datasets)) == 0:
                logger.debug('Ignoring dataset not in list of datasets to consider: "%s"' % ds_name)
                self._ignored_datasets.append(ds_name)
                return None

            # Add a new instance of the dataset from a new file
            ds_info.add_instance(hdf_obj, ds_obj)

            if len(ds_info.inp_shape_names) > 0:
                self.datasets_info[ds_name] = ds_info
            else:
                logger.warning( "Ignoring dataset which is not an array: %s" % ds_name)

        hdf_obj.visititems(_check_dataset)

    def _sort_file_lists(self):
        # Sort filenames, sounding ids and indexes based on sounding ids
        sort_lists = zip(self.file_sounding_ids, self.file_indexes, self.src_filenames)
        sort_lists.sort(key=itemgetter(0)) # Sort by sounding id
        self.file_sounding_ids = map(itemgetter(0), sort_lists)
        self.file_indexes = map(itemgetter(1), sort_lists)
        self.src_filenames = map(itemgetter(2), sort_lists)

        # Make a list if getter found only one item and made it just a single string
        if not hasattr(self.src_filenames, "__iter__"):
            self.file_sounding_ids = [ self.file_sounding_ids ]
            self.file_indexes = [ self.file_indexes ]
            self.src_filenames = [ self.src_filenames ]

    def analyze_files(self):
        # Loop over a copy of the source filenames as we will remove from the original list any that
        # have no matching sounding ids
        check_files = copy(self.src_filenames)

        for file_idx, curr_file in enumerate(check_files):
            with closing(InputContainerChooser(curr_file, single_id_dim=not self.multi_source_types)) as hdf_obj:
                num_ids = self.process_file_ids(hdf_obj)

                if num_ids > 0:
                    self.process_dataset_info(hdf_obj)
                else:
                    self.src_filenames.remove(curr_file)

            write_progress(file_idx, len(check_files))
    
        self._sort_file_lists()

def create_output_dataset(out_hdf_obj, dataset_info, splice_size=None, collapse_id_dim=False):
    """Duplicates a dataset from the input file into the output hdf object as it exists
    except for its dimensions"""

    logger.debug("Creating new output dataset: %s" % dataset_info.out_name)

    dst_shape, max_shape, dst_shape_names = dataset_info.output_dataset_shape(splice_size, collapse_id_dim) 

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
        out_dataset_obj = out_group_obj.create_dataset(dst_name, shape=dst_shape, dtype=dataset_info.out_type, maxshape=max_shape, compression="gzip", compression_opts=2, fillvalue=dataset_fill)
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
    if dst_shape_names:
        out_dataset_obj.attrs["Shape"] = numpy.array(["_".join(dst_shape_names) + "_Array"]) 

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

def create_dest_file_datasets(out_hdf_obj, source_info, copy_all=False, multi_source_types=False):

    # Loop over dataset names/shapes
    # Create datasets and copy any non id_shape datasets
    num_datasets = len(source_info.datasets_info.keys())
    out_dataset_objs = {}
    for curr_index, (curr_dataset_name, curr_dataset_info) in enumerate(source_info.datasets_info.items()):
        write_progress(curr_index-1, num_datasets)

        id_dimension = curr_dataset_info.inp_id_dims

        if out_hdf_obj.get(curr_dataset_info.out_name, None) != None:
            logger.debug("Destination dataset %s for %s already exists in output file, ignoring duplicate" % (curr_dataset_info.out_name, curr_dataset_name))
            del source_info.datasets_info[curr_dataset_name]

        elif len(id_dimension) > 0: 
            # Create a dataset that does have an id dimension
            out_dataset_objs[curr_dataset_name] = create_output_dataset(out_hdf_obj, curr_dataset_info, source_info.num_soundings, collapse_id_dim=multi_source_types)
        
        elif copy_all:
            logger.debug("Non sounding id dataset, copying from first file: %s" % curr_dataset_name)
                
            # Copy dataset with no sounding dimension directly from first file
            out_dataset_obj = create_output_dataset(out_hdf_obj, curr_dataset_info, collapse_id_dim=multi_source_types)
           
            # Search for first file that has the desired dataset. Can't just take
            # from first in cases when splicing dissimilar files
            for inp_fn in source_info.src_filenames:
                with closing(h5py.File(inp_fn, 'r')) as hdf_obj:
                    if hdf_obj.get(curr_dataset_name, None):
                        # Do not try and copy empty datasets
                        if hdf_obj[curr_dataset_name].len() > 0:
                            out_dataset_obj[:] = hdf_obj[curr_dataset_name][:]
                        source_info.datasets_info.pop(curr_dataset_name)
                        break
        else:
            del source_info.datasets_info[curr_dataset_name]
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

def get_datasets_for_file(curr_file, curr_snd_indexes, inp_datasets_info, output_index, multi_source_types=False):
    # Ignore current file if no soundings
    if len(curr_snd_indexes) == 0:
        return     

    with closing(InputContainerChooser(curr_file, single_id_dim=not multi_source_types)) as curr_hdf_obj:
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
                        out_dataset_idx, out_shape = output_indexes_shape(curr_ds_obj, curr_dataset_info, output_slice, collapse_id_dim=multi_source_types)

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

def copy_multiple_datasets(out_dataset_objs, source_info, multi_source_types=False):
    """Given an input and output hdf object, copies the specified soundings ids for the specified
    dataset names along the specified shape dimention to the output file"""

    # Loop over dataset names/shapes
    output_index = 0
    for file_index, (curr_file, curr_snd_indexes) in enumerate(zip(source_info.src_filenames, source_info.file_indexes)):
        logger.debug( 'Copying %d sounding(s) from %s (%d / %d)' % (len(curr_snd_indexes), os.path.basename(curr_file), file_index+1, len(source_info.src_filenames)) )

        write_progress(file_index-1, len(source_info.src_filenames))

        for curr_dataset_name, out_dataset_idx, stored_data in get_datasets_for_file(curr_file, curr_snd_indexes, source_info.datasets_info, output_index, multi_source_types):
            out_dataset_objs[curr_dataset_name].__setitem__(tuple(out_dataset_idx), stored_data)

        write_progress(file_index, len(source_info.src_filenames))

        # Increment where to start in the index into output sounding indexes
        # Only increment if not using multiple source types because when using different
        # inputs to a single output we will be pulling the same sounding ids (hopefully)
        if not multi_source_types:
            output_index += len(curr_snd_indexes) 

def determine_multi_source(check_files):

    # Looks through input hdf files and see if we are using multiple source types, ie different types of files combined into one
    multi_source_types = False

    file_id_names = []
    for curr_file in check_files:
        with closing(InputContainerChooser(curr_file)) as hdf_obj:
            file_id_names.append( hdf_obj.get_id_dim_names() )

    # Lets make sure the inputs are either all different types
    # or all the same, otherwise if some are the same,
    # bad things could happen
    id_match_counts = numpy.array([ len(filter(lambda x: x == idn, file_id_names)) for idn in file_id_names ])

    if(all(id_match_counts == len(id_match_counts))):
        logger.debug("Not using multiple source types")
    elif(len(id_match_counts) > 1 and all(id_match_counts >= 1) and all(id_match_counts < len(id_match_counts))):
        multi_source_types = True
        logger.debug("Using multiple source types")
    else:
        raise Exception("Error determining if using multiple source types with id match counts: %s" % id_match_counts)

    return multi_source_types

def splice_files(source_files, output_filename, sounding_ids, splice_all=False, desired_datasets_list=None, rename_mapping=False, multi_source_types=None):

    logger.info("Splicing into: %s" % output_filename)

    # Determine if input files are of different types
    if multi_source_types == None:
        logger.info("Determining if using multiple source types")
        multi_source_types = determine_multi_source(source_files)

    # Match sounding ids to files and indexes using L1B files
    source_info = SourceInformation(source_files, sounding_ids, desired_datasets_list, rename_mapping, multi_source_types)
    with Timer() as t:
        logger.info("Analyzing source files")
        source_info.analyze_files()

    logger.info("Analyzed %d files, found %d sounding ids and %d unique datasets in %d files taking %.03f seconds" % (len(source_files), len(source_info.found_sounding_ids), len(source_info.datasets_info), len(source_info.src_filenames), t.interval))

    # Copy over all datasets in files
    with closing(h5py.File(output_filename, 'w')) as dest_obj:

        # Create datasets in the output file so we have somewhere to dump the input
        with Timer() as t:
            logger.info("Creating output file datasets")
            out_dataset_objs = create_dest_file_datasets(dest_obj, source_info, splice_all, multi_source_types)
            dest_obj.flush()

        logger.info("Datasets creation took %.03f seconds" % (t.interval))

        # Now copy desired datasets from source file to destination files
        with Timer() as t:
            logger.info("Copying data to output")
            copy_multiple_datasets(out_dataset_objs, source_info, multi_source_types=multi_source_types)

        logger.info("Datasets copying took %.03f seconds" % (t.interval))

def load_source_files(filenames, input_list_file=None):
    source_list = []
    if len(filenames) > 0:
        source_list += filenames

    if input_list_file != None:
        source_list += [ fn.strip() for fn in open(input_list_file).readlines() ]
    
    return source_list

def standalone_main():
    parser = ArgumentParser(description="Splices together input HDF5 products for a given set of sounding ids")

    parser.add_argument( "filenames", metavar="FILE", nargs='*',
                         help="files to splice, may be left blank if using the -i --input-files-list option" )

    parser.add_argument( "-i", "--input-files-list", dest="input_files_list",
                         metavar="FILE",
                         help="text file with input filenames to splice")

    parser.add_argument( "-s", "--sounding-id-file", dest="sounding_id_file",
                         metavar="FILE",
                         help="file containing list of soundings for destination file")

    parser.add_argument( "-o", "--output-file", dest="output_filename",
                         metavar="FILE", default=DEFAULT_OUTPUT_FILENAME,
                         help="output filename of splice data, default: %s" % DEFAULT_OUTPUT_FILENAME)

    parser.add_argument( "-d", "--datasets-list-file", dest="datasets_list_file",
                         metavar="FILE",
                         help="file containing list of only datasets to consider for copying. If rename_mapping is enabled then the names are matched on their destination dataset name")

    parser.add_argument( "-r", "--rename-mapping", dest="rename_mapping",
                         action="store_true",
                         default=False,
                        help="rename datasets into output file according to internal mapping table as they would appear in the L2Agg PGE")

    parser.add_argument( "--agg-names-filter", dest="agg_names_filter",
                         action="store_true",
                         default=False,
                         help="include only dataset names that would appear in the L2Agg PGE. Its only makes sense to use this option with --rename_mapping")

    parser.add_argument( "--splice-all", dest="splice_all",
                         action="store_true",
                         default=False,
                         help="splice all datasets, including those which do not have a sounding dimension. Note that datasets without an explicit handler and no sounding dimension are simply copied from the first file.")

    parser.add_argument( "--multiple-file-types", dest="multi_source_types", action="store_true", default=None,
                         help="indicates that multiple file type sources are being spliced. Speeds up multiple source type determination stage by being specified." )

    parser.add_argument( "--single-file-type", dest="multi_source_types", action="store_false", default=None,
                         help="indicates that a single type of file is being spliced. Speeds up multiple source type determination stage by being specified." )

    parser.add_argument( "-v", "--verbose", dest="verbose",
                         action="store_true",
                         default=False,
                         help="enable verbose informational reporting")

    # Parse command line arguments
    args = parser.parse_args()

    if len(args.filenames) == 0 and args.input_files_list == None:
        parser.error('input list file must be specified')

    # Set up logging
    if args.verbose:
        # Include HDF5 errors in output
        h5py._errors.unsilence_errors()

        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    console = logging.StreamHandler(stream=sys.stdout)
    logger.addHandler(console)

    source_files = load_source_files(args.filenames, args.input_files_list)

    if args.sounding_id_file != None:
        sounding_ids = [ sid.strip() for sid in open(args.sounding_id_file).readlines() ]
    else:
        logger.info("No sounding ids file supplied, aggregating all ids from all files")
        sounding_ids = None

    copy_datasets_list = None
    if args.datasets_list_file != None:
        copy_datasets_list = [ ds.strip() for ds in open(args.datasets_list_file).readlines() ]

    if args.agg_names_filter:
        if copy_datasets_list == None:
            copy_datasets_list = aggregator_dataset_dest_names
        else:
            copy_datasets_list += aggregator_dataset_dest_names

    splice_files(source_files, args.output_filename, sounding_ids, splice_all=args.splice_all, desired_datasets_list=copy_datasets_list, rename_mapping=args.rename_mapping, multi_source_types=args.multi_source_types)

if __name__ == "__main__":
    standalone_main()
