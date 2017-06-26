#!/usr/bin/env python2.3

import sys, os, re
from mechanize import Browser
from time import strftime, strptime
import ClientForm

####
# Check command line arguments

if (len(sys.argv) < 7):
    print "usage:\n\t", os.path.basename(sys.argv[0]), "<location_id> <lat_beg> <lat_end> <lon_beg> <lon_end> <date_beg> <date_end>\n"
    print "location_id is any string without whitespace"
    print "latitude as DDD[N|S], ex: 90N or 90S"
    print "longitude as DDD[W|E],  ex: 180W or 180E"
    print "time as YYYYMMDD, ex: 20080801"
    print ""
    print "example:\n\t", os.path.basename(sys.argv[0]), "JPL 30N 35N 240E 245E 20080101 20080131"
    sys.exit(1)

(location_id, lat_beg, lat_end, lon_beg, lon_end, date_beg, date_end) = sys.argv[1:8]

if location_id.find(' ') >= 0:
    print 'location_id can not have a space in it'
    sys.exit(1)

if not re.match('\d+\.?(\d+)?[NS]', lat_beg):
    print 'lat_beg: "%s" not in proper format: DDD[N|S]' % lat_beg

if not re.match('\d+\.?(\d+)?[NS]', lat_end):
    print 'lat_end "%s" not in proper format: DDD[N|S]' % lat_end

if not re.match('\d+\.?(\d+)?[EW]', lon_beg):
    print 'lon_beg "%s" not in proper format: DDD[E|W]' % lon_beg

if not re.match('\d+\.?(\d+)?[EW]', lon_end):
    print 'lon_end "%s" not in proper format: DDD[E|W]' % lon_end

if not re.match('\d{8}', date_beg):
    print 'date_beg "%s" no in proper format: YYYYMMDD' % date_beg

if not re.match('\d{8}', date_end):
    print 'date_end "%s" no in proper format: YYYYMMDD' % date_end

print 'arguments =', location_id, lat_beg, lat_end, lon_beg, lon_end, date_beg, date_end

###
# Parse date in proper format
try:
    time_beg_struct = strptime(date_beg, '%Y%m%d')
except ValueError, e_obj:
    print 'error parsing date_beg "%s" : %s' % (date_beg, e_obj)
    sys.exit(1)
try:
    time_end_struct = strptime(date_end, '%Y%m%d')
except ValueError, e_obj:
    print 'error parsing date_end "%s" : %s' % (date_end, e_obj)
    sys.exit(1)

year_beg = strftime('%Y', time_beg_struct)
year_end = strftime('%Y', time_end_struct)

mon_beg = strftime('%b', time_beg_struct)
mon_end = strftime('%b', time_end_struct)

day_beg = '%d' % int(strftime('%d', time_beg_struct))
day_end = '%d' % int(strftime('%d', time_end_struct))

print 'date_beg =', year_beg, mon_beg, day_beg
print 'date_end =', year_end, mon_end, day_end

#sys.exit(1)

# eg..  NCEP_JPL_30N_35N_240E_245E_20080201_20080229_6hourly_AT.nc
# location_id, lat_beg, lat_end, lon_beg, lon_end, date_beg, date_end, set
out_file_format = 'NCEP_%s_%s_%s_%s_%s_%s_%s_6hourly_%s.nc'

browse = Browser()

# Masquerade as a real person and not a robot
browse.set_handle_robots(False)

ncep_set_urls = {
    'AT' : 'http://www.cdc.noaa.gov/cgi-bin/db_search/DBSearch.pl?Dataset=NCEP+Reanalysis+Pressure+Level&Dataset=NCEP+Reanalysis+Daily+Averages+&Pressure+Level&Variable=Air+Temperature',
    'GH' : 'http://www.cdc.noaa.gov/cgi-bin/db_search/DBSearch.pl?Dataset=NCEP+Reanalysis+Pressure+Level&Dataset=NCEP+Reanalysis+Daily+Averages+Pressure+Level&Variable=Geopotential+height',
    'SH' : 'http://www.cdc.noaa.gov/cgi-bin/db_search/DBSearch.pl?Dataset=NCEP+Reanalysis+Pressure+Level&Dataset=NCEP+Reanalysis+Daily+Averages+&Pressure+Level&Variable=Specific+humidity',
    'TP' : 'http://www.cdc.noaa.gov/cgi-bin/db_search/DBSearch.pl?Dataset=NCEP+Reanalysis+Tropopause+Level&Dataset=NCEP+Reanalysis+Daily+Averages+Tropopause+Level&Variable=Pressure',
    }

for (curr_set_name, curr_set_url) in ncep_set_urls.iteritems():
    print 'Loading page for set: %s' % curr_set_name
    print 'URL: %s' % curr_set_url
    
    browse.open(curr_set_url)
    browse.follow_link(url_regex=re.compile(r"DB_statistic\=Individual"))
    
    browse.select_form(nr=1)

    try:
        level_obj = browse.form.find_control("level")
        has_levels = True
    except ClientForm.ControlNotFoundError:
        has_levels = False
    
    if has_levels:
        all_levels = []
        for lev_info in level_obj.get_items():
            all_levels.append(lev_info.name)
        print 'Levels = ', all_levels
        browse["level"]  = all_levels

    browse["lon-begin"] = lon_beg
    browse["lon-end"]   = lon_end
    browse["lat-begin"] = lat_beg
    browse["lat-end"]   = lat_end
    
    browse["year_begin"] = [year_beg]
    browse["mon_begin"]  = [mon_beg]
    browse["day_begin"]  = [day_beg]
    browse["hour_begin"] = ["00 Z"]

    browse["year_end"]   = [year_end]
    browse["mon_end"]    = [mon_end]
    browse["day_end"]    = [day_end]
    browse["hour_end"]   = ["18 Z"]

    browse["output"] = ["file"]

    print 'Submitting form for set: %s' % curr_set_name 
    browse.submit()

    out_filename = out_file_format % (location_id, lat_beg, lat_end, lon_beg, lon_end, date_beg, date_end, curr_set_name)

    #print 'Links:'
    #for curr_link in browse.links():
    #    print curr_link
    
    request = browse.follow_link(url_regex=re.compile(r"ftp\:\/\/"))
    print 'Downloading from url: %s' % browse.geturl()
    if os.path.exists(out_filename):
        print("%s already exists, not grabbing" % out_filename)
    else:
        print 'Saving to: %s' % out_filename
        file_obj = file(out_filename, "wb")
        while 1:
            data = request.read(1024)
            if not data: break
            sys.stdout.write('.')
            file_obj.write(data)
        print ''
        file_obj.close()

