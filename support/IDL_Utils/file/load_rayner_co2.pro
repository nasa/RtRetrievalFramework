pro load_rayner_co2, month, pressures, longitudes, latitudes, julian_times, file_co2

    co2_files_glob_fmt='("/data2/groups/algorithm/mcduffie/rayner_4hb/cpickett.lmdz.histrac.an2000.m",I02,".nc")'

    latitudes=[90.0000, 87.5000, 85.0000, 82.5000, 80.0000, 77.5000, 75.0000, 72.5000, 70.0000, 67.5000, 65.0000, 62.5000, 60.0000, 57.5000, 55.0000, 52.5000, 50.0000, 47.5000, 45.0000, 42.5000, 40.0000, 37.5000, 35.0000, 32.5000, 30.0000, 27.5000, 25.0000, 22.5000, 20.0000, 17.5000, 15.0000, 12.5000, 10.0000, 7.50000, 5.00000, 2.50000, 0.00000, -2.50000, -5.00000, -7.50000, -10.0000, -12.5000, -15.0000, -17.5000, -20.0000, -22.5000, -25.0000, -27.5000, -30.0000, -32.5000, -35.0000, -37.5000, -40.0000, -42.5000, -45.0000, -47.5000, -50.0000, -52.5000, -55.0000, -57.5000, -60.0000, -62.5000, -65.0000, -67.5000, -70.0000, -72.5000, -75.0000, -77.5000, -80.0000, -82.5000, -85.0000, -87.5000, -90.00]

    longitudes=[-180.000, -176.250, -172.500, -168.750, -165.000, -161.250, -157.500, -153.750, -150.000, -146.250, -142.500, -138.750, -135.000, -131.250, -127.500, -123.750, -120.000, -116.250, -112.500, -108.750, -105.000, -101.250, -97.5000, -93.7500, -90.0000, -86.2500, -82.5000, -78.7500, -75.0000, -71.2500, -67.5000, -63.7500, -60.0000, -56.2500, -52.5000, -48.7500, -45.0000, -41.2500, -37.5000, -33.7500, -30.0000, -26.2500, -22.5000, -18.7500, -15.0000, -11.2500, -7.50000, -3.75000,  0.00000,  3.75000,  7.50000,  11.2500,  15.0000,  18.7500,  22.5000,  26.2500,  30.0000,  33.7500,  37.5000,  41.2500,  45.0000,  48.7500,  52.5000,  56.2500,  60.0000,  63.7500,  67.5000,  71.2500,  75.0000,  78.7500,  82.5000,  86.2500,  90.0000,  93.7500,  97.5000,  101.250,  105.000,  108.750,  112.500,  116.250,  120.000,  123.750,  127.500,  131.250,  135.000,  138.750,  142.500,  146.250,  150.000,  153.750,  157.500,  161.250,  165.000,  168.750,  172.500,  176.250]

    ;; Load filenames
    search_glob=string(month, format=co2_files_glob_fmt)
    input_file=file_search(search_glob)
    
    if input_file[0] eq '' then $
      Message, 'Could not find file: ', search_glob
       
    print, 'Loading: ', input_file[0]
    ncid = ncdf_open(input_file[0])
    
    ;; Pressure column same for all files
    print, 'Reading pressure axis'
    ncdf_varget, ncid, 'presnivs', pressures
    
    ;; Get starting time for time axis
    ncdf_attget, ncid, 'time_counter', 'units', time_units
    time_units = string(time_units) ;; byte array to string
    
    if strpos(time_units, 'seconds since ') ne 0 then $
      Message, 'Units for time does not have expected format'
    
    beg_year=fix(strmid(time_units, strlen('seconds since '), 4))
    beg_month=fix(strmid(time_units, strlen('seconds since ')+5,2))
    beg_day=fix(strmid(time_units, strlen('seconds since ')+8, 2))
    
    ;; Convert times to julian 
    print, 'Reading time axis'
    ncdf_varget, ncid, 'time_counter', file_times
    julian_times=timegen(n_elements(file_times), units="Days", start=julday(beg_month, beg_day, beg_year, 00, 00, 00), seconds=file_times)
    
    print, 'Reading co2 values'
    ncdf_varget, ncid, 'mthly_3rd', file_co2
    
    ncdf_close, ncid

    help, pressures, longitudes, latitudes
    help, julian_times, file_co2

end
