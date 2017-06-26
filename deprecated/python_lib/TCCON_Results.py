import numpy
from datetime import datetime, timedelta

class TCCON_File(object):
    def __init__(self):
        raise Exception("Abstract base class only")

    def read_column_data(self, all_columns, data_row, header_skip):
        # Search for columns we can not parse as numbers and seperate
        # from those we can not
        self.columns = []
        addl_col_names = []

        data_col_indexes = []
        addl_col_indexes = []
        for col_idx, data_val in enumerate(data_row.strip().split()):
            try:
                float_val = float(data_val.strip())
                self.columns.append(all_columns[col_idx])
                data_col_indexes.append(col_idx)
            except ValueError:
                addl_col_names.append(all_columns[col_idx])
                addl_col_indexes.append(col_idx)

        # Parse the parts of the file that can be used to create a big matrix
        self.data = numpy.loadtxt(self.filename, skiprows=header_skip, usecols=data_col_indexes)

        # Go through and assign values in data to individual attributes
        for data_idx, data_name in enumerate(self.columns):
            setattr(self, data_name, self.data[:, data_idx])

        # Add additional non-float types as attributes
        addl_values = numpy.loadtxt(self.filename, skiprows=header_skip, usecols=addl_col_indexes, dtype=str)
        # Make sure we have a 2d array even if just 1 column
        if len(addl_values.shape) == 1:
            addl_values = addl_values.reshape((addl_values.shape[0], 1))

        for data_idx, data_name in enumerate(addl_col_names):
            setattr(self, data_name, addl_values[:, data_idx])

        # Process times from runlog into a datetime objects
        self.process_datetimes()

    def process_datetimes(self):
        self.datetimes = []
        if hasattr(self, "year") and hasattr(self, "day") and hasattr(self, "hour"):
            for year, doy, tday in zip(self.year, self.day, self.hour):
                hour = int(tday)
                if hour >= 24:
                    hour = 0
                    tday -= 24
                    doy += 1

                minute = int((tday - hour) * 60.0)
                if minute >= 60:
                    minute -= 60
                    hour += 1

                second = int(((tday - hour) * 60.0 - minute) * 60.0)
                frac_sec = ((tday - hour) * 60.0 - minute) * 60.0 - second
                micro_second = int(1e6 * frac_sec)

                dt = datetime(int(year), 1, 1, hour, minute, second, micro_second) + timedelta(days=int(doy-1))
                self.datetimes.append(dt)

class Database_Results(TCCON_File):
    """Reads data from a table file downloaded from http://tccon.ipac.caltech.edu/"""
    
    def __init__(self, filename=None):
        if filename != None:
            self.read(filename)

    def read(self, filename):
        # Assign filename here in case a different file is read later on
        self.filename = filename

        # Gather column names
        line_count = 0
        with open(filename, "r") as data_obj:
            # Find column line
            line = data_obj.readline()
            while(len(line) > 0):
                line_count += 1
                if line.find("|") == 0:
                    column_line = line
                    break
                line = data_obj.readline()
        
            # Strip off newline then | at beginning and end, split by | and remove extra
            # whitespace in the column names
            all_columns = [ c.strip() for c in column_line.strip().strip("|").split("|") ]

            # Find first data row
            while(len(line) > 0 and line.find("|") == 0):
                line_count += 1
                line = data_obj.readline()

        self.read_column_data(all_columns, line, line_count-1)

class GGG_Output(TCCON_File):
    """Reads a file written by GGG"""

    def __init__(self, filename=None):
        if filename != None:
            self.read(filename)

    def read(self, filename):
        # Assign filename here in case a different file is read later on
        self.filename = filename

        with open(filename, "r") as f_obj:
            header_line = f_obj.readline().strip()
            try:
                (num_skip, num_cols, num_rows) = [ int(val) for val in header_line.strip().split() ][:3]
            except:
                raise ValueError("Could not parse num_skip, num_rows, num_cols from header line: " + header_line)

            header_len = num_skip - 1 
            while (header_len > 0):
                header_line = f_obj.readline().strip()
                header_len -= 1

            all_columns = header_line.split()

            first_data_line = f_obj.readline().strip()

        self.read_column_data(all_columns, first_data_line, num_skip)

    def write_fake_l2_file(self, output_filename, xco2=None):
        sounding_ids = numpy.array([ int(s[-3:]) for s in self.spectrum ])

        if xco2 == None:
            if hasattr(self, "xco2"):
                xco2 = self.xco2
            else:
                raise Exception("GGG file %s does not have an xco2 column" % self.filename)

        import h5py
        h_out = h5py.File(output_filename, "w")

        shape_attr = numpy.array(["Retrieval_Array"])

        sid_ds = h_out.create_dataset("/RetrievalHeader/sounding_id_reference", data=sounding_ids)
        sid_ds.attrs['Shape'] = shape_attr

        xco2_ds = h_out.create_dataset("/RetrievalResults/xco2", data=xco2)
        xco2_ds.attrs['Shape'] = shape_attr
        
        time_strings = numpy.array([ [dt.strftime("%Y-%m-%dT%H:%M:%S") for x in range(3)] for dt in self.datetimes ])
        time_string_ds = h_out.create_dataset("/FtsRunLog/time_string", data=time_strings)
        time_string_ds.attrs['Shape'] = numpy.array(["Retrieval_Dim_1"])
        time_string_ds.attrs['Type'] = numpy.array(["VarLenStr"])

        solar_zeniths = numpy.array([ [s for x in range(3)] for s in self.asza ])
        solar_zenith_ds = h_out.create_dataset("/FtsRunLog/solar_zenith", data=solar_zeniths)
        solar_zenith_ds.attrs['Shape'] = numpy.array(["Retrieval_Dim_1"])

        h_out.close()

class VAV_File(GGG_Output):
    """Reads a VAV file, and allows for creation of a fake L2 file for plotting purposes"""

    def write_fake_l2_file(self, output_filename):
        xco2 = self.co2 / (self.o2 / 0.2095)
        GGG_Output.write_fake_l2_file(self, output_filename, xco2)

