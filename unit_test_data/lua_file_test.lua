test_data_dir = os.getenv("abs_top_srcdir") .. "/unit_test_data/"
sid_string = "20090725020316"

l1b_fname = test_data_dir .. "l1b.h5"
l1b_hdf = HdfFile(l1b_fname)
sid = AcosSoundingId(l1b_hdf, sid_string, AcosSoundingId.P_SOUNDING)


