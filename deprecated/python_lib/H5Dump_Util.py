import os
import sys
from OCO_MathUtil import *
from types import ListType

h5dump_bin = "h5dump"

class H5Dump_Util:

    def __init__(self, filename, debug=False):

        self.filename = filename
        self.debug = debug
        
        self.dims = []

        if not os.path.exists(filename):
            raise IOError('Could not find HDF file: %s' % filename)

    def get_all_dataset_names(self):
        h5dump_cmd = "%s -n %s" % (h5dump_bin, self.filename)

        h5dump_fo = os.popen(h5dump_cmd)
        in_contents = False
        dataset_names = []
        for h5dump_line in h5dump_fo.readlines():
            h5dump_line = h5dump_line.strip()

            if self.debug:
                print 'Parsing line: %s' % h5dump_line
            
            if h5dump_line.find('FILE_CONTENTS') == 0:
                in_contents = True
            elif h5dump_line.find('}') == 0:
                in_contents = False
                break
            elif in_contents == True:
                (ds_type, ds_name) = h5dump_line.split()
                if ds_type == 'dataset':
                    dataset_names.append(ds_name)
            
        h5dump_fo.close()

        return dataset_names

    def get_dims(self, dataset):
        h5dump_cmd = "%s -A -d %s %s" % (h5dump_bin, dataset, self.filename)

        if self.debug:
            print 'Running command: %s' % h5dump_cmd

        dims = None
        h5dump_fo = os.popen(h5dump_cmd)
        for h5dump_line in h5dump_fo.readlines():
            h5dump_line = h5dump_line.strip()

            if h5dump_line.find('DATASPACE') >= 0:
                if self.debug:
                    print 'Found DATASPACE line:\n-> %s' % h5dump_line

                dim_beg = h5dump_line.find('{')+1
                dim_end = h5dump_line.rfind('}')

                dim_sects = (h5dump_line[dim_beg+1:dim_end]).split('/')
                
                if len(dim_sects) < 1:
                    raise IOError('Could not find dimensions section for dataset %s in file %s' % (dataset, self.filename))
                
                dims = [ int(dval) for dval in dim_sects[0].replace('(','').replace(')','').split(',') ]
                break

        if dims == None:
            raise IOError('Could not parse dimensions for dataset %s in file %s' % (dataset, self.filename))
        
        h5dump_fo.close() 

        return dims

    def get_dataset_values(self, dataset, start=None, stride=None, count=None, ignore_set_error=False, datatype=float):

        dims = self.get_dims(dataset)

        opt_command_hash = {'start':'-s', 'stride':'-S', 'count':'-c'}
        additional_opts = []
        for opt_name in ['start', 'stride', 'count']:
            opt_cmd = opt_command_hash[opt_name]
            opt_vals = eval(opt_name)
            if opt_vals != None and len(opt_vals) != len(dims):
                raise ValueError('Length of %s keyword: %d must match dimensions of dataset: %d' % (opt_name, len(opt_vals), len(dims)))
            elif opt_vals != None:
                opt_str = '%s %s' % ( opt_cmd, ','.join([ str(val) for val in opt_vals ]) )
                additional_opts.append(opt_str)
                                   
        h5dump_cmd = "%s -w 15 -d %s %s %s" % (h5dump_bin, dataset, ' '.join(additional_opts), self.filename)

        if self.debug:
            print 'Running command: %s' % h5dump_cmd

        data_count = dims
        dataset_values = None
        extra_set_idx = 0
            
        h5dump_fo = os.popen(h5dump_cmd)
        in_data_sect = False
        for h5dump_line in h5dump_fo.readlines():
            h5dump_line = h5dump_line.strip()

            if self.debug:
                print 'Parsing line: %s' % h5dump_line

            if h5dump_line.find('DATA {') == 0:
                in_data_sect = True
                dataset_values = numpy.zeros(data_count, dtype=datatype)
            elif in_data_sect and h5dump_line.find('}') == 0:
                in_data_sect = False
                break
            elif h5dump_line.find('COUNT (') == 0:
                count_beg = h5dump_line.find('(')
                count_end = h5dump_line.rfind(')')
                data_count = [ int(cstr) for cstr in h5dump_line[count_beg+1:count_end].split(',') ]
                print data_count
                for cnt_idx in range(len(data_count)):
                    if data_count[cnt_idx] > 1:
                        extra_set_idx = cnt_idx
                        break
            elif in_data_sect:
                (idx_str, data_str) = [part.strip() for part in h5dump_line.split(':')]

                idx_beg = idx_str.find('(')
                idx_end = idx_str.rfind(')')
                
                indexes = [ int(istr) for istr in idx_str[idx_beg+1:idx_end].split(',') ]
                data_strs = data_str.rstrip(',').split(',')

                if dataset_values == None:
                    raise ValueError('dataset_values not yet initialized')

                if len(indexes) != len(data_count):
                    raise ValueError('Number of dimensions for data value do not match declared dimensions')
                    
                for c_idx in range(len(data_count)):
                    if indexes[c_idx] > (data_count[c_idx] - 1):
                        raise ValueError('Index %s outside that of allocated data ( %s )' % (idx_str, ', '.join(data_count)))

                str_idx = 0
                index_locations = [ int(idx) for idx in idx_str[idx_beg+1:idx_end].split(',') ]
                for curr_str in data_strs:
                    index_locations[extra_set_idx] += str_idx
                    set_command = 'dataset_values[%s] = %s' % (','.join([str(idx) for idx in index_locations]), curr_str)
                    try:
                        exec(set_command)
                    except:
                        print >>sys.stderr, 'Error setting data values:'
                        print >>sys.stderr, 'Index String:\n->%s' % idx_str
                        if not ignore_set_error:
                            print >>sys.stderr, 'Data line:\n->%s' % h5dump_line
                            print >>sys.stderr, 'Data String:\n->%s' % data_str
                            print >>sys.stderr, 'Set Command:\n->%s' % set_command
                            raise ValueError('Error setting data values with set command: %s' % set_command)
                    str_idx += 1

        return dataset_values
