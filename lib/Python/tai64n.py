""" Module to read and convert TAI64(n) hex strings as used by DJBDNS.

    decode_tai64n("hex_string")
        decodes a hex string representing a TAI64(n) external date and returns
        a datetime.datetime object at UTC.
    utc2tai(datetime.datetime.now())
        returns a datetime.datetime object representing TAI from a datetime.datetime
        object at UTC. 
"""
from __future__ import division
from past.utils import old_div


from datetime import datetime, timedelta
from operator import itemgetter


class TAI64DecodeError(Exception):
    """ TAI64 decode error.
    """
    pass


def __conversion_table():
    """ returns [datetime, value] ordered reverse by date
        where value == seconds between TAI and UTC
    
    Example:
    >>> __conversion_table()[-1]
    (datetime.datetime(1972, 1, 1, 0, 0), 10.0)
    """
    # update this table as new values become known
    # source: ftp://maia.usno.navy.mil/ser7/tai-utc.dat
    conversion_table = [(datetime(1972, 0o1,  1), 10.0),
                        (datetime(1972, 0o7,  1), 11.0),
                        (datetime(1973, 0o1,  1), 12.0),
                        (datetime(1974, 0o1,  1), 13.0),
                        (datetime(1975, 0o1,  1), 14.0),
                        (datetime(1976, 0o1,  1), 15.0),
                        (datetime(1977, 0o1,  1), 16.0),
                        (datetime(1978, 0o1,  1), 17.0),
                        (datetime(1979, 0o1,  1), 18.0),
                        (datetime(1980, 0o1,  1), 19.0),
                        (datetime(1981, 0o7,  1), 20.0),
                        (datetime(1982, 0o7,  1), 21.0),
                        (datetime(1983, 0o7,  1), 22.0),
                        (datetime(1985, 0o7,  1), 23.0),
                        (datetime(1988, 0o1,  1), 24.0),
                        (datetime(1990, 0o1,  1), 25.0),
                        (datetime(1991, 0o1,  1), 26.0),
                        (datetime(1992, 0o7,  1), 27.0),
                        (datetime(1993, 0o7,  1), 28.0),
                        (datetime(1994, 0o7,  1), 29.0),
                        (datetime(1996, 0o1,  1), 30.0),
                        (datetime(1997, 0o7,  1), 31.0),
                        (datetime(1999, 0o1,  1), 32.0),
                        (datetime(2006, 0o1,  1), 33.0),
                        (datetime(2009, 0o1,  1), 34.0),
                        # add new values here
                       ]
    conversion_table.sort(key=itemgetter(0), reverse=True)
    return conversion_table

def __tai_seconds(date, table=__conversion_table()):
    """ returns seconds of TAI-offset from UTC at date given.
        Works only on dates later than 01.01.1972.
    
    Example:
        >>> __tai_seconds(datetime(1992, 6, 2, 8, 7, 9))
        26.0
        >>> __tai_seconds(datetime(1971, 6, 2, 8, 7, 9))
        False
    """
    for x in table:
        if date > x[0]:
            return x[1]
    return False


def tai2utc(date):
    """ converts datetime.datetime TAI to datetime.datetime UTC.
        Works only on dates later than 01.01.1972.
    
    Example
        >>> tai2utc(datetime(1992, 6, 2, 8, 7, 9))
        datetime.datetime(1992, 6, 2, 8, 6, 43)
    """
    seconds = __tai_seconds(date)
    return seconds and (date - timedelta(0, seconds))


def utc2tai(date):
    """ converts datetime.datetime UTC to datetime.datetime TAI.
        Works only on dates later than 01.01.1972.

    Example
        >>> utc2tai(datetime(1992, 6, 2, 8, 6, 43))
        datetime.datetime(1992, 6, 2, 8, 7, 9)
    """
    seconds = __tai_seconds(date)
    return seconds and (date + timedelta(0, seconds))


def decode_tai64n(hexstring, basedate=datetime(1970, 1, 1)):
    """ returns a datetime.datetime (UTC) object from a tai64(n) string.
        Works only on dates after 01.01.1972.
        
        Args
            hexstring: string containig TAI64(n) value
            basedate: DO NOT USE (speeds up repeated use of this function)
        
        Returns
            a datetime.datetime object at UTC
        
        Example
            (from http://cr.yp.to.mirror.dogmap.org/libtai/tai64.html)
            >>> decode_tai64n("400000002a2b2c2d")
            datetime.datetime(1992, 6, 2, 8, 6, 43)
            >>> decode_tai64n("400000004a32392b2aa21dac")
            datetime.datetime(2009, 6, 12, 11, 16, 25, 715267)
    """
    try:
        nano_hex = (len(hexstring) == 24) and hexstring[16:24]
    except IndexError:
        pass
    try:
        tai_int = int(hexstring[0:16], 16)
        nano_int = nano_hex and int(nano_hex, 16)
    except:
        raise TAI64DecodeError("'%s' not a valid hex value." % hexstring)
    # we decode only dates later than 01.01.1972
    seconds = tai_int - 4611686018427387904
    if seconds < 0:
        raise TAI64DecodeError("I won't decode gone millenia "
                               "(i.e. nothing prior to 01.01.1972).")
    return tai2utc(basedate + timedelta(0, seconds, old_div(nano_int,1000)))



