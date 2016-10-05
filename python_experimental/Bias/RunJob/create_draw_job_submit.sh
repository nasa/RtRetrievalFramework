#!/bin/bash
numsounding=`cat sounding_list.txt | wc -w`
numdraw=700
endindex=`expr $numsounding \* $numdraw - 1`
qsub -j oe -o log/create_draw_job.log -t 0-${endindex} ./create_draw_job.sh -q long
