# Test ILS type permutations

begin KEYWORDS
   configuration_sub_dirs    = <configuration_sub_dirs>/<ils_type_name>
end KEYWORDS

# LIST defining to turn LSI off and on
begin LIST
   name = ils_type_name
   begin VALUES
      sinc
      gauss
      table
   end VALUES
end LIST
