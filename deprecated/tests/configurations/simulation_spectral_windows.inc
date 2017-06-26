# Should be included within a run directory config file.
# Sets up keywords for observation mode spectral windows

begin MAP
    name = window_map_1
    begin VALUES
       gosat     ABO2_Obs
       oco       ''
       oco_tconn ''
    end VALUES
end MAP

begin MAP
    name = window_map_2
    begin VALUES
       gosat     WCO2_Obs
       oco       ''
       oco_tconn ''
    end VALUES
end MAP

begin MAP
    name = window_map_3
    begin VALUES
       gosat     SCO2_Obs
       oco       ''
       oco_tconn ''
    end VALUES
end MAP

begin KEYWORDS
    spectral_windows_1 = [window_map_1:instrument_name]
    spectral_windows_2 = [window_map_2:instrument_name]
    spectral_windows_3 = [window_map_3:instrument_name]
end KEYWORDS
