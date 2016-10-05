begin FILE
   source_filename = <case_psurf_filename>

   #######
   # PSURF

   begin MODIFY_DATA
      columns = PSURF
      modify  = 'abs( <original> + random.gauss(0,200) )'
   end MODIFY_DATA

end FILE
