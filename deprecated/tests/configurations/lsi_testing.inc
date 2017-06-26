# Test permutations of lsi on, off

begin KEYWORDS
   configuration_sub_dirs    = <configuration_sub_dirs>/[lsi_name:use_lsi]
end KEYWORDS

# LIST defining to turn LSI off and on
begin LIST
   name = use_lsi
   begin VALUES
      True
      False
   end VALUES
end LIST

# MAP translating true/false for LSI into strings for directory name
begin MAP
   name = lsi_name
   begin VALUES
      True  lsi_on
      False lsi_off
   end VALUES
end MAP
