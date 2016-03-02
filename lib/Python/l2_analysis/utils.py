from __future__ import division
from past.utils import old_div
import time
from datetime import datetime, timedelta

def id2time(sounding_ids):
    times = []
    for snd_id in sounding_ids:
        try:
            snd_str = "%s" % snd_id
            dt = datetime.strptime(snd_str[:14], "%Y%m%d%H%M%S")
        except TypeError as exc:
            raise TypeError("Could not parse sounding id string: %s" % snd_id)

        if len(snd_str) > 14:
            frac_sec = snd_str[14]
            if frac_sec == "0":
                pass
            elif frac_sec == "3" or frac_sec == "4":
                dt += timedelta(seconds=old_div(1.0,3.0))
            elif frac_sec == "7":
                dt += timedelta(seconds=old_div(2.0,3.0))
            else:
                raise Exception("Unknown fractional seconds integer: %s in sounding id: %s" % (snd_id, frac_sec))
        times.append(dt)
    return times

def time_string_to_dt(time_strings):
    return [ datetime.strptime(ts[:19], "%Y-%m-%dT%H:%M:%S") for ts in time_strings ] 
