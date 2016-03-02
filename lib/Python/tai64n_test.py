from __future__ import absolute_import
from nose.tools import *
from .tai64n import *
from datetime import datetime, timedelta

def test_tai2utc():
    assert tai2utc(datetime(1992, 6, 2, 8, 7, 9)) == datetime(1992, 6, 2, 8, 6, 43)

def test_utc2tai():
    assert utc2tai(datetime(1992, 6, 2, 8, 6, 43)) == datetime(1992, 6, 2, 8, 7, 9)

def test_decode_tai64n():
    assert decode_tai64n("400000002a2b2c2d") == datetime(1992, 6, 2, 8, 6, 43)
    assert decode_tai64n("400000004a32392b2aa21dac") == datetime(2009, 6, 12, 11, 16, 25, 715267)
