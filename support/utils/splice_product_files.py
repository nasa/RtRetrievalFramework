#!/usr/bin/env python
# This utility takes a file containing a list of L1B/ECMWF/Cloud filenames as well as a file
# containing a list of sounding ids. The program will then extract the data for those sounding
# ids from the various files and splice them into a new L1B and and additional files.

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import zip
from builtins import map
from builtins import range
from builtins import object

import os
import re
import sys
import math
from copy import copy
from argparse import ArgumentParser
from operator import itemgetter
from collections import Counter, Iterable
from contextlib import closing
import logging
import traceback
import time

import multiprocessing
from tempfile import mkstemp

# Should be progressbar2
from progressbar import ProgressBar, Percentage, SimpleProgress, Bar

# For terminal control
from blessings import Terminal

import numpy
import h5py

import full_physics.acos_file as acos_file
from full_physics.fill_value import FILL_VALUE
from full_physics import log_util

from full_physics.splice_tool_mapping import aggregator_dataset_mapping, aggregator_dataset_dest_names

DEFAULT_OUTPUT_FILENAME = "spliced_output.h5"

# Base progress bar widgets
PROGRESS_WIDGETS = [ Bar(), ' ', Percentage(), ' (', SimpleProgress(), ')' ]

term = Terminal()

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

class Timer(object):    
    def __enter__(self):
        self.start = time.process_time()
        return self

    def __exit__(self, *args):
        self.end = time.process_time()
        self.interval = self.end - self.start

class ProgressWriter(object):
    "Writes progress bar to a specific location on the screen, but only if terminal supports it, other wise no progress bar is printed"

    def __init__(self, location=(), fd=sys.stdout):
        """
        Input: location - tuple of ints (x, y), the position
                          of the bar in the terminal
        """

        self.location = location
        self.fd = fd

    def write(self, string):
        # Only write progress is terminal supports it and is not being piped
        if term.does_styling:
            if self.location != None and len(self.location) > 0:
                with term.location(*self.location):
                    print(string, file=self.fd)
            else:
                self.fd.write(string)

    def flush(self):
        self.fd.flush()

class DatasetInformation(object):

    def __init__(self, dataset_name, rename_mapping=False, collapse_id_dim=False):
        # Strip leading and trailing slashes just in case
        self.inp_name = dataset_name.strip("/")

        base_name = self.inp_name.split("/")[-1]
        if rename_mapping and base_name in aggregator_dataset_mapping:
            self.out_name = aggregator_dataset_mapping[base_name]["name"]
            self._override_shape_names = aggregator_dataset_mapping[base_name]["shape"].replace("_Array", "").split("_")
        else:
            self.out_name = self.inp_name
            self._override_shape_names = None

        self.collapse_id_dim = collapse_id_dim

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

        # Eliminate duplicate shape names
        for name, num_occurance in list(Counter(self.inp_shape_names).items()):
            if num_occurance > 1:
                for name_count in range(1, num_occurance+1):
                    self.inp_shape_names[self.inp_shape_names.index(name)] = "%s%d" % (name, name_count)

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

    @property
    def out_shape_names(self):
        if self._override_shape_names is not None:
            out_shape_names = self._override_shape_names
        else:
            out_shape_names = list(self.inp_shape_names)

        # Concatenate input shape names for id dims into one name
        # Go through and collapse self.inp_id_dims dim(s) into a single one of
        # sized according to self.inp_id_dims
        if self.collapse_id_dim and self._override_shape_names is None and out_shape_names is not None and len(self.inp_id_dims) > 0:
            out_id_name = ''.join([ out_shape_names[didx] for didx in self.inp_id_dims ])
            out_shape_names[self.inp_id_dims[0]] = out_id_name

            # Remove any other splice dims,
            # Iterate backwards since we are
            for s_dim in self.inp_id_dims[1:][::-1]:
                out_shape_names.pop(s_dim)

        return out_shape_names

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
            raise ValueError('No shape names are defined for dataset: %s' % (dataset_name))
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


        return tuple(data_dim_indexes)
 
    def output_dataset_shape(self, splice_size=None):
        # Make a copy of dataset shape size and names as they will be modified
        dst_shape = self.out_shape 

        if self.inp_id_dims != None and len(self.inp_id_dims) > 0 and splice_size != None:
            # Set first splice dim to splice size
            dst_shape[self.inp_id_dims[0]] = splice_size

        if self.collapse_id_dim: 
            # Remove any other splice dims,
            # Iterate backwards since we are
            # popping from dst_shape
            for s_dim in self.inp_id_dims[1:][::-1]:
                dst_shape.pop(s_dim)

        # Make sure any zero sized datasets are indicated as unlimited
        # by identifying them as None in the maxshape argument
        max_shape = [None if s == 0 else s for s in dst_shape]

        return dst_shape, max_shape

class SourceInformation(object):

    def __init__(self, src_filenames, sounding_ids=None, desired_datasets=None, rename_mapping=False, multi_source_types=False, logger=logging.getLogger(), progress=ProgressBar()):

        # Filenames to splice from, will be modified to remove any which have no sounding ids present
        # Ensure it is a list so we can modify it as we go
        self.src_filenames = list(src_filenames)

        if sounding_ids != None:
            self.inp_sounding_ids = set(numpy.array(sounding_ids, dtype=numpy.int64))
        else:
            self.inp_sounding_ids = None

        # All sounding ids found from source files that match input list
        self.found_sounding_ids = set() # Will be a list after analyze is run
        self.found_ids_to_index = {}

        # Flags
        self.desired_datasets = desired_datasets
        self.rename_mapping = rename_mapping
        self.multi_source_types = multi_source_types

        # Progress bar
        self.progress = copy(progress)
        self.progress.widgets = [ 'Analyzing sources ' ] + self.progress.widgets

        # Logger to use for class
        self.logger = logger

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
        curr_sounding_ids = numpy.ravel(hdf_obj.get_sounding_ids())

        # Convert curr_sounding_ids to a set in place, because if we use all sounding ids
        # we don't want an unsorted set used which can mess up the ordering of sounding ids
        if self.inp_sounding_ids == None:
            matching_ids = curr_sounding_ids
        else:
            matching_ids = sorted(set(curr_sounding_ids).intersection(self.inp_sounding_ids))

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
                ds_info = DatasetInformation(ds_name, rename_mapping=self.rename_mapping, collapse_id_dim=self.multi_source_types)

            # Search list of datasets based upon the destination dataset name, which is the same as
            # the input name when rename mapping is not used
            if self.desired_datasets != None and len([ds for ds in self.desired_datasets if re.search(ds_info.out_name, ds)]) == 0:
                self.logger.debug('Ignoring dataset not in list of datasets to consider: "%s"' % ds_name)
                self._ignored_datasets.append(ds_name)
                return None

            # Add a new instance of the dataset from a new file
            ds_info.add_instance(hdf_obj, ds_obj)

            if len(ds_info.inp_shape_names) > 0:
                self.datasets_info[ds_name] = ds_info
            else:
                self.logger.warning( "Ignoring dataset which is not an array: %s" % ds_name)

        hdf_obj.visititems(_check_dataset)

    def _sort_file_lists(self):
        # Sort filenames, sounding ids and indexes based on sounding ids
        sort_lists = list(zip(self.file_sounding_ids, self.file_indexes, self.src_filenames))
        sort_lists.sort(key=itemgetter(0)) # Sort by sounding id
        self.file_sounding_ids = list(map(itemgetter(0), sort_lists))
        self.file_indexes = list(map(itemgetter(1), sort_lists))
        self.src_filenames = list(map(itemgetter(2), sort_lists))

        # Make a list if getter found only one item and made it just a single string
        if not hasattr(self.src_filenames, "__iter__"):
            self.file_sounding_ids = [ self.file_sounding_ids ]
            self.file_indexes = [ self.file_indexes ]
            self.src_filenames = [ self.src_filenames ]

        # Make sure list of destination ids is in a sorted order
        self.found_sounding_ids = sorted(self.found_sounding_ids)

        # Make a mapping of sounding ids to index
        self.found_ids_to_index = { sid: index[0] for (index, sid) in numpy.ndenumerate(self.found_sounding_ids) }

    def analyze_files(self):
        # Loop over a copy of the source filenames as we will remove from the original list any that
        # have no matching sounding ids
        check_files = copy(self.src_filenames)

        self.progress.start(len(check_files))
        for file_idx, curr_file in enumerate(check_files):
            with closing(InputContainerChooser(curr_file, single_id_dim=not self.multi_source_types)) as hdf_obj:
                num_ids = self.process_file_ids(hdf_obj)

                if num_ids > 0:
                    self.process_dataset_info(hdf_obj)
                else:
                    self.src_filenames.remove(curr_file)

            self.progress.update(file_idx + 1)

        self.progress.finish()
        self.progress.fd.write('\n')
        
        self._sort_file_lists()

class ProductSplicer(object):

    def __init__(self, source_files, output_filename, sounding_ids, splice_all=False, desired_datasets_list=None, rename_mapping=False, multi_source_types=None, logger=logging.getLogger(), progress=ProgressBar()):
        
        self.splice_all = splice_all
        self.multi_source_types = multi_source_types

        self.progress = copy(progress)
        self.progress.widgets = ['Copying data '] + self.progress.widgets

        self.logger = logger

        self.source_info = SourceInformation(source_files, sounding_ids, desired_datasets_list, rename_mapping, multi_source_types, logger, progress)
        self.dest_obj = h5py.File(output_filename, 'w')
        self.out_dataset_objs = {}

    def __exit__(self):
        "Close output file on loss of context"
        self.dest_obj.close()

    def analyze(self):

        self.source_info.analyze_files()
    
    def create_output_dataset(self, dataset_info, splice_size=None):
        """Duplicates a dataset from the input file into the output hdf object as it exists
        except for its dimensions"""

        self.logger.debug("Creating new output dataset: %s" % dataset_info.out_name)

        dst_shape, max_shape = dataset_info.output_dataset_shape(splice_size) 

        # Split name into two and create a group if needed
        # then create the desire dataset or load existing group
        ds_name_clean = dataset_info.out_name.lstrip('/').rstrip('/')
        if ds_name_clean.find("/") > 0:
            dst_group, dst_name = dataset_info.out_name.lstrip('/').rstrip('/').split('/', 1)
            out_group_obj = self.dest_obj.require_group(dst_group) 
        else:
            dst_group = ""
            dst_name = ds_name_clean
            out_group_obj = self.dest_obj

        # Fill new dataset with the correct fill value based on type
        if dataset_info.out_type != numpy.object and dataset_info.out_type.type != numpy.string_:
            fill_type = dataset_info.out_type.type

            if fill_type in FILL_VALUE:
                dataset_fill = FILL_VALUE[fill_type]
            else:
                self.logger.warning("Could not find specific fill value for dataset: %s of type %s" % (dst_name, fill_type))
                dataset_fill = None
        else:
            # Use default fill for string types
            fill_type = None
            dataset_fill = None

        self.logger.debug( "Creating new dataset: %s/%s sized: %s with fill type: %s and value: %s" % (dst_group, dst_name, dst_shape, fill_type, dataset_fill) )
        try:
            out_dataset_obj = out_group_obj.create_dataset(dst_name, shape=dst_shape, dtype=dataset_info.out_type, maxshape=max_shape, compression="gzip", compression_opts=2, fillvalue=dataset_fill, chunks=True)
        except RuntimeError as exc:
            raise RuntimeError("Error creating dataset %s/%s: %s" % (dst_group/dst_name, exc))

        # Now create copied attributes from original dataset
        # Just copy from first for now, leave code to do multiple if needed
        for curr_file in dataset_info.inp_filenames[0:1]:
            with closing(h5py.File(curr_file, 'r')) as curr_hdf_obj:
                curr_dataset_obj = curr_hdf_obj[dataset_info.inp_name]
                for attr_name, attr_value in list(curr_dataset_obj.attrs.items()):
                    # Skip if copied already from a file, assuming all files
                    # have same attributes for now, we just will get
                    # all uniquely named ones
                    if attr_name in list(out_dataset_obj.attrs.keys()):
                        continue
            
                    # If the dtype of the attribute is an object dtype, then assume
                    # its a variable length string
                    if hasattr(attr_value, "dtype") and attr_value.dtype.kind == "O":
                        self.logger.debug('Copying variable length string attribute: "%s" with value: "%s"' % (attr_name, attr_value[0]))
                        vlen_dt = h5py.special_dtype(vlen=str)
                        out_dataset_obj.attrs.create(attr_name, attr_value, dtype=vlen_dt)
                    else:
                        self.logger.debug('Copying attribute: "%s" with value: "%s"' % (attr_name, attr_value))
                        out_dataset_obj.attrs.create(attr_name, attr_value)

        # Add extra information for dataset, overwrite an existing shape, because we may have
        # reshaped the data. This must be a bytes and not string object
        if dataset_info.out_shape_names is not None:
            out_dataset_obj.attrs["Shape"] = numpy.array([("_".join(dataset_info.out_shape_names) + "_Array").encode('utf-8')]) 

        return out_dataset_obj

    def create_datasets(self):

        # Loop over dataset names/shapes
        # Create datasets and copy any non id_shape datasets
        num_datasets = len(list(self.source_info.datasets_info.keys()))
        for curr_index, (curr_dataset_name, curr_dataset_info) in enumerate(list(self.source_info.datasets_info.items())):
            id_dimension = curr_dataset_info.inp_id_dims

            if self.dest_obj.get(curr_dataset_info.out_name, None) != None:
                self.logger.debug("Destination dataset %s for %s already exists in output file, ignoring duplicate" % (curr_dataset_info.out_name, curr_dataset_name))
                del self.source_info.datasets_info[curr_dataset_name]

            elif len(id_dimension) > 0: 
                # Create a dataset that does have an id dimension
                self.out_dataset_objs[curr_dataset_name] = self.create_output_dataset(curr_dataset_info, self.source_info.num_soundings)
            
            elif self.splice_all:
                self.logger.debug("Non sounding id dataset, copying from first file: %s" % curr_dataset_name)
                    
                # Copy dataset with no sounding dimension directly from first file
                out_dataset_obj = self.create_output_dataset(curr_dataset_info)
               
                # Search for first file that has the desired dataset. Can't just take
                # from first in cases when splicing dissimilar files
                for inp_fn in self.source_info.src_filenames:
                    with closing(h5py.File(inp_fn, 'r')) as hdf_obj:
                        if hdf_obj.get(curr_dataset_name, None):
                            # Do not try and copy empty datasets
                            if hdf_obj[curr_dataset_name].len() > 0:
                                out_dataset_obj[:] = hdf_obj[curr_dataset_name][:]
                            self.source_info.datasets_info.pop(curr_dataset_name)
                            break
            else:
                del self.source_info.datasets_info[curr_dataset_name]
                self.logger.debug("Ignoring non sounding id dataset: %s" % curr_dataset_name)

        self.dest_obj.flush()

    def output_indexes_shape(self, in_dataset_obj, dataset_info, out_data_idx):
        # Create output dataset indexes, removing any concatenated dimensions 
        splice_dim = dataset_info.inp_id_dims 
        out_dataset_idx = [ slice(dim_size) for dim_size in in_dataset_obj.shape ]
        out_dataset_idx[splice_dim[0]] = out_data_idx

        if self.multi_source_types:
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
            elif isinstance(idx, Iterable):
                out_shape.append(len(idx))
            elif type(idx) is int:
                out_shape.append(1)
            else:
                raise ValueError("Unknown type %s in indexes: %s" % (type(idx), out_dataset_idx))

        return out_dataset_idx, out_shape

    def get_dataset_slice(self, in_dataset_obj, dataset_info, in_data_idx, out_shape, inp_filename=""):
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
            raise IOError("Can not read dataset %s from file %s: %s" % (dataset_info.inp_name, inp_filename, exc))

        # Set sliced data into output dataset
        if numpy.product(in_data.shape) > numpy.product(out_shape):
            self.logger.warning("Dataset %s requires destructive resizing" % (dataset_info.out_name))
            self.logger.debug("At indexes %s resizing source data of shape %s to %s." % (in_data_idx, in_data.shape, out_shape)) 
            stored_data = numpy.resize(in_data, out_shape)
        else:
            stored_data = in_data.reshape(out_shape)

        return stored_data

    def get_datasets_for_file(self, curr_file, curr_snd_indexes, output_slice):
        # Ignore current file if no soundings
        if len(curr_snd_indexes) == 0:
            return     

        with closing(InputContainerChooser(curr_file, single_id_dim=not self.multi_source_types)) as curr_hdf_obj:
            # Copy each dataset for current file
            for curr_dataset_name, curr_dataset_info in list(self.source_info.datasets_info.items()):

                curr_ds_obj = curr_hdf_obj.get(curr_dataset_name, None)
                if curr_ds_obj != None:
                    # Only actually perform the copy if the object does not contain an empty dimension
                    # otherwise we leave the created output dataset empty
                    if not 0 in curr_ds_obj.shape:
                        # Copy dataset's relevant contents into output
                        try:
                            out_dataset_idx, out_shape = self.output_indexes_shape(curr_ds_obj, curr_dataset_info, output_slice)

                            stored_data = self.get_dataset_slice(curr_ds_obj, curr_dataset_info, curr_snd_indexes, out_shape, curr_hdf_obj.filename)
                            yield curr_dataset_name, out_dataset_idx, stored_data
          
                        except (ValueError, IOError) as e:
                            self.logger.error('Error copying dataset: "%s", dataset may be corrupt: %s' % (curr_dataset_name, e))
                            self.logger.debug('Exception: %s' % traceback.format_exc())
                    else:
                        self.logger.debug('Dataset: "%s" has an empty dimension in its shape: %s, not actually copying any values from: %s' % (curr_dataset_name, curr_ds_obj.shape, curr_file))

    def copy_datasets(self):
        """Given an input and output hdf object, copies the specified soundings ids for the specified
        dataset names along the specified shape dimention to the output file"""

        self.progress.start(len(self.source_info.src_filenames))

        # For each file copy all of its datasets to the destination file
        for file_index, (curr_file, curr_snd_ids, curr_snd_indexes) in enumerate(zip(self.source_info.src_filenames, self.source_info.file_sounding_ids, self.source_info.file_indexes)):
            # Output as simpler message on progress when progress bar output is disabled
            if not term.does_styling:
                self.logger.info('Copying %d sounding(s) from %s (%d / %d)' % (len(curr_snd_indexes), os.path.basename(curr_file), file_index+1, len(self.source_info.src_filenames)))

            # Figure out where into the output file to copy source data for this file
            output_indexes = numpy.array([ self.source_info.found_ids_to_index[sid] for sid in curr_snd_ids ])
            index_count = numpy.sum(output_indexes[1:] - output_indexes[:-1]) + 1

            if index_count != output_indexes.shape[0]:
                # Slices are not contiguous so copy each index individually, slower
                output_slices = []
                for indiv_in_index, indiv_out_index in zip(curr_snd_indexes, output_indexes):
                    output_slices.append( ([indiv_in_index], slice(indiv_out_index, indiv_out_index + 1)) )
            else:
                # When output indexes are contiguous, faster
                output_slices = [ (curr_snd_indexes, slice(output_indexes[0], output_indexes[-1] + 1)) ]

            for curr_in_indexes, curr_out_slice in output_slices:
                for curr_dataset_name, out_dataset_idx, stored_data in self.get_datasets_for_file(curr_file, curr_in_indexes, curr_out_slice):
                    self.out_dataset_objs[curr_dataset_name].__setitem__(tuple(out_dataset_idx), stored_data)

            self.progress.update(file_index)

        self.progress.finish()
        self.progress.fd.write('\n')

        # Flush our copies to disk
        self.dest_obj.flush()

    def splice_files(self):

        self.logger.info("Splicing into: %s" % self.dest_obj.filename)

        # Match sounding ids to files and indexes using L1B files
        with Timer() as t:
            self.logger.info("Analyzing source files")
            self.analyze()

        self.logger.info("Found %d sounding ids and %d unique datasets in %d files taking %.03f seconds" % (len(self.source_info.found_sounding_ids), len(self.source_info.datasets_info), len(self.source_info.src_filenames), t.interval))

        # Create datasets in the output file so we have somewhere to dump the input
        with Timer() as t:
            self.logger.info("Creating output file datasets")
            self.create_datasets()

        self.logger.info("Datasets creation took %.03f seconds" % (t.interval))

        # Now copy desired datasets from source file to destination files
        with Timer() as t:
            self.logger.info("Copying data to output")
            self.copy_datasets()

        self.logger.info("Datasets copying took %.03f seconds" % (t.interval))
    
def determine_multi_source(check_files, logger):

    # Looks through input hdf files and see if we are using multiple source types, ie different types of files combined into one
    multi_source_types = False

    file_id_names = []
    for curr_file in check_files:
        with closing(InputContainerChooser(curr_file)) as hdf_obj:
            file_id_names.append( hdf_obj.get_id_dim_names() )

    # Lets make sure the inputs are either all different types
    # or all the same, otherwise if some are the same,
    # bad things could happen
    id_match_counts = numpy.array([ len([x for x in file_id_names if x == idn]) for idn in file_id_names ])

    if(all(id_match_counts == len(id_match_counts))):
        logger.debug("Not using multiple source types")
    elif(len(id_match_counts) > 1 and all(id_match_counts >= 1) and all(id_match_counts < len(id_match_counts))):
        multi_source_types = True
        logger.debug("Using multiple source types")
    else:
        raise Exception("Error determining if using multiple source types with id match counts: %s" % id_match_counts)

    return multi_source_types


def splice_worker(worker_idx, source_files, output_fn, sounding_ids, splice_all, desired_datasets_list, rename_mapping, multi_source_types, log_file):
    # Put worker names into logger and progress bars
    progress = ProgressBar(widgets=PROGRESS_WIDGETS, fd=ProgressWriter((0, worker_idx)))
    progress.widgets.append(' - Worker #%02d' % (worker_idx + 1))

    logger = logging.getLogger("Worker%02d" % (worker_idx + 1))
    logger.propagate = False
    if log_file != None:
        logger.addHandler(log_file)
    else:
        # Use a null handler to avoid no handlers could be found message
        class NullHandler(logging.Handler):
            def emit(self, record):
                pass
        logger.addHandler(NullHandler())

    worker_splicer = ProductSplicer(source_files, output_fn, sounding_ids, splice_all, desired_datasets_list, rename_mapping, multi_source_types, logger, progress)
    worker_splicer.splice_files()

def process_files(source_files, output_filename, sounding_ids, splice_all=False, desired_datasets_list=None, rename_mapping=False, multi_source_types=None, workers=1, temp_dir=None, main_logger=logging.getLogger(), log_file=None):

    # Sort source files to reduce chance of out of order sounding ids when processing in parallel
    source_files.sort()

    # Determine if input files are of different types
    if multi_source_types == None:
        main_logger.info("Determining if using multiple source types")
        multi_source_types = determine_multi_source(source_files, main_logger)

    if workers > 1 and len(source_files) >= workers:
        # Clear screen to make cursor indexing easier
        sys.stdout.write(term.clear)
    
        # Make space for progress bars
        sys.stdout.write(term.move_down * workers)

        # Split source files into equal groups, one group for each worker
        group_size = int(math.ceil(len(source_files) / workers))
        source_groups = [ source_files[pos:pos + group_size] for pos in range(0, len(source_files), group_size) ]

        main_logger.info("Splicing into %d temporary files in parallel" % len(source_groups))

        processes = []
        final_files = []
        temp_files = []

        # Launch workers then wait for their completion
        for worker_idx in range(len(source_groups)):
            # Save into a temporary file
            temp_obj, temp_fn = mkstemp(suffix=".h5", dir=temp_dir)

            final_files.append(temp_fn)
            temp_files.append(temp_fn)

            proc = multiprocessing.Process(target=splice_worker, args=(worker_idx, source_groups[worker_idx], temp_fn, sounding_ids, splice_all, desired_datasets_list, rename_mapping, multi_source_types, log_file))
            processes.append(proc)

        # Start processes
        for proc in processes:
            proc.start()

        # Wait till all processes end
        for proc in processes:
            proc.join()

        # Clear screen for final splicing if just used parallel mode
        sys.stdout.write(term.clear)
        final_writer = ProgressWriter((0, 0))

    else:
        # No parallelization, source files ae the final files we deal with
        final_files = source_files
        temp_files = []
        final_writer = ProgressWriter()

    progress = ProgressBar(widgets=PROGRESS_WIDGETS, fd=final_writer)

    final_splicer = ProductSplicer(final_files, output_filename, sounding_ids, splice_all, desired_datasets_list, rename_mapping, multi_source_types, main_logger, progress)
    final_splicer.splice_files()

    for fn in temp_files:
        os.unlink(fn)

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

    parser.add_argument( "-w", "--workers", dest="workers", type=int, default=1,
                         help="Number of workers to use when parallelizing splicing" )

    parser.add_argument( "--temp", dest="temp_dir", default=os.curdir,
                         help="Directory where temporary files are saved when number of parallel workers is greater than 1" )

    parser.add_argument( "-l", "--log_file", dest="log_file", 
                         help="Save verbose information to log file" )

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
        parser.error("Input list file must be specified")

    # Set up logging
    if args.verbose:
        # Include HDF5 errors in output
        h5py._errors.unsilence_errors()

        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    main_logger = log_util.init_logging(log_level=log_level, format="%(message)s")
        
    # Initialize logging
    if args.log_file:
        log_file = log_util.open_log_file(args.log_file, logger=main_logger)
        log_file.setFormatter( logging.Formatter("%(asctime)s: %(name)8s - %(levelname)7s - %(message)s") )
    else:
        log_file = None

    source_files = load_source_files(args.filenames, args.input_files_list)

    if args.sounding_id_file != None:
        sounding_ids = [ sid.strip() for sid in open(args.sounding_id_file).readlines() ]
    else:
        main_logger.debug("No sounding ids file supplied, aggregating all ids from all files")
        sounding_ids = None

    copy_datasets_list = None
    if args.datasets_list_file != None:
        copy_datasets_list = [ ds.strip() for ds in open(args.datasets_list_file).readlines() ]

    if args.agg_names_filter:
        if copy_datasets_list == None:
            copy_datasets_list = aggregator_dataset_dest_names
        else:
            copy_datasets_list += aggregator_dataset_dest_names

    process_files(source_files, args.output_filename, sounding_ids, splice_all=args.splice_all, desired_datasets_list=copy_datasets_list, rename_mapping=args.rename_mapping, multi_source_types=args.multi_source_types, workers=args.workers, temp_dir=args.temp_dir, main_logger=main_logger, log_file=log_file)

if __name__ == "__main__":
    standalone_main()
