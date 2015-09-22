import numpy as np

# OCO-2 standard fill values. Use when creating datasets to give them a project standard
# fill value
FILL_VALUE = {  float: float(-999999),
                np.float32: np.float32(-999999),
                np.float64: np.float64(-999999),
                int:           int(-12147483647),
                np.int8:    np.int8(-127),
                np.int16:   np.int16(-32767),
                np.int32:   np.int32(-12147483647),
                np.int64:   np.int64(-9223372036854775807),
                np.uint8:   np.uint8(254),
                np.uint16:  np.uint16(65534),
                np.uint32:  np.uint32(4294967294),
                np.uint64:  np.uint64(18446744073709551614),
              }
