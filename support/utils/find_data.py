#!/usr/bin/env python
from __future__ import division
from builtins import str
from past.utils import old_div
import sqlite3
from pyproj import Geod

# Distance in km
def distance(lat1, lon1, lat2, lon2):
    g = Geod(ellps='WGS84')
    return old_div(g.inv(lon1, lat1, lon2, lat2)[2],1000)

# alter table l1b add column julian_time real;
# update l1b set julian_time = julianday(time);
# create index time_index on l1b(julian_time);
conn = sqlite3.connect("l1b.db")
conn.create_function("distance", 4, distance)

c = conn.cursor()
date1 = '2010-05-15'
date2 = '2010-06-15'
lat = 33.7443
lon = 118.3491
dist = 500
file = {}
sid_file = open("sid.txt", "w")
for row in c.execute("select sounding_id,file,ecmwf from l1b where julian_time >= julianday(?) and julian_time <= julianday(?) and distance(latitude, longitude, ?, ?) < ?", (date1, date2, lat, lon, dist)):
    sid_file.write(str(row[0]) + "\n")
    file[row[1]] = row[2]

input = open("input.txt", "w")
for f,e in file.items():
    input.write(str(f) + " " + str(e) + "\n")


