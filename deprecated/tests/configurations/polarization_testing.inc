# Test permutations of polarization on, off

begin KEYWORDS
   configuration_sub_dirs    = <configuration_sub_dirs>/[pol_name:use_polarization]
end KEYWORDS

# LIST defining to turn polarization off and on
begin LIST
   name = use_polarization
   begin VALUES
      True
      False
   end VALUES
end LIST

# MAP translating true/false for polarization into strings for directory name
begin MAP
   name = pol_name
   begin VALUES
      True  pol_on
      False pol_off
   end VALUES
end MAP
