#!/usr/bin/env python

# Add some fields to spliced L2 files which maps the L1B indexes into where they are used in the file
# Requires two argments, the L1B filename and the spliced file to write into

import os
import sys
import h5py
from bisect import bisect

from argparse import ArgumentParser

import numpy 

from full_physics.fill_value import FILL_VALUE

IN_L1B_SND_ID_DS = "/SoundingGeometry/sounding_id"
RET_SND_IDS_DS = "/RetrievalHeader/sounding_id"

# Maps L1B indexes to retrieval indexes
OUT_L1B_SND_ID_DS = "/L1bScSoundingReference/sounding_id_l1b"
RET_INDEX_DS = "/L1bScSoundingReference/retrieval_index"

# Maps retrieval indexes to L1B indexes
OUT_FRAME_INDEX_DS = "/RetrievalHeader/frame_index"
OUT_SOUNDING_INDEX_DS = "/RetrievalHeader/sounding_index"

# Datasets that may be missing from the output that should be added as empty datasets
ADD_EMPTY_DATASETS = {
        ("/RetrievalResults/wind_speed", numpy.float64),
        ("/RetrievalResults/wind_speed_apriori", numpy.float64),
        ("/RetrievalResults/wind_speed_uncert", numpy.float64),
        ("/RetrievalResults/fluorescence_at_reference", numpy.float64),
        ("/RetrievalResults/fluorescence_at_reference_apriori", numpy.float64),
        ("/RetrievalResults/fluorescence_at_reference_uncert", numpy.float64),
        ("/RetrievalResults/fluorescence_slope", numpy.float64),
        ("/RetrievalResults/fluorescence_slope_apriori", numpy.float64),
        ("/RetrievalResults/fluorescence_slope_uncert", numpy.float64),
        }

parser = ArgumentParser(description='Add fields into a spliced retrieval output file mapping L1B soundings to Retrieval indexes')
parser.add_argument('l1b_filename', help='Input L1B filename to generate mapping from')
parser.add_argument('splice_filename', help='Spliced retrieval file to place mapping inside')

args = parser.parse_args()

l1b_file = h5py.File(args.l1b_filename, "r")
splice_file = h5py.File(args.splice_filename, "r+")

l1b_sounding_ids = l1b_file.get(IN_L1B_SND_ID_DS, None)
if l1b_sounding_ids == None:
    raise Exception("Could not find sounding id data set from l1b file: %s" % args.l1b_filename)

ret_sounding_ids = splice_file.get(RET_SND_IDS_DS, None)
if ret_sounding_ids == None:
    raise Exception("Could not find retrieval sounding ids data set from spliced file: %s" % args.splice_filename)

def copy_attrs(in_ds, out_ds):
    for k, v in in_ds.attrs.items():
        out_ds.attrs[k] = numpy.array([v[0]])

# Copy L1B sounding ids into destination file
if splice_file.get(OUT_L1B_SND_ID_DS, None) == None:
    out_snd_ids_ds = splice_file.create_dataset(OUT_L1B_SND_ID_DS, data=l1b_sounding_ids)
    copy_attrs(l1b_sounding_ids, out_snd_ids_ds)
else:
    print >>sys.stderr, "Dataset %s already exists in %s" % (OUT_L1B_SND_ID_DS, args.splice_filename)

# Add index datasets that index from L1B file indexes to L2 file indexes and vice versa
if splice_file.get(RET_INDEX_DS, None) == None:
    retrieval_index = numpy.zeros(l1b_sounding_ids.shape, dtype=int)

    frame_index = numpy.zeros(ret_sounding_ids.shape, dtype=int)
    sounding_index = numpy.zeros(ret_sounding_ids.shape, dtype=int)

    for ret_idx, ret_id in enumerate(ret_sounding_ids):
        where_id = numpy.where(l1b_sounding_ids[:, :] == ret_id)

        if len(where_id) > 0:
            frame_idx = where_id[0][0]
            snd_idx = where_id[1][0]

            if int(l1b_sounding_ids[frame_idx, snd_idx]) != int(ret_id):
                raise Exception("Found L1B sounding id: %d which is not the same as retrieval sounding id used to search: %d at frame index: %d" % (l1b_sounding_ids[frame_idx, snd_idx], ret_id, frame_idx))

            retrieval_index[frame_idx, snd_idx] = ret_idx+1

            frame_index[ret_idx] = frame_idx + 1
            sounding_index[ret_idx] = snd_idx + 1

    out_retrieval_index_ds = splice_file.create_dataset(RET_INDEX_DS, data=retrieval_index)
    copy_attrs(l1b_sounding_ids, out_retrieval_index_ds)

    out_frame_index_ds = splice_file.create_dataset(OUT_FRAME_INDEX_DS, data=frame_index)
    copy_attrs(ret_sounding_ids, out_frame_index_ds)

    out_sounding_index_ds = splice_file.create_dataset(OUT_SOUNDING_INDEX_DS, data=sounding_index)
    copy_attrs(ret_sounding_ids, out_sounding_index_ds)
else:
    print >>sys.stderr, "Dataset %s already exists in %s" % (RET_INDEX_DS, args.splice_filename)

# Add datasets that may be missing due to them being specific to, say a ground type, but the aggregator never
# encountered any soundings that had the dataset
for ds_name, ds_type in ADD_EMPTY_DATASETS:
    if splice_file.get(ds_name, None) == None:
        new_ds = splice_file.create_dataset(ds_name, data=numpy.empty(ret_sounding_ids.shape, dtype=ds_type))
        new_ds[:] = FILL_VALUE.get(ds_type, FILL_VALUE[float])

l1b_file.close()
splice_file.close()
