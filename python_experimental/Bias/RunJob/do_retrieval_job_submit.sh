#!/bin/bash
numsounding=`cat sounding_list.txt | wc -w`
numdraw=700
endindex=`expr $numsounding \* $numdraw - 1`
qsub -j oe -o log/retrieval_initial.log -t 0-${endindex} ./do_retrieval_job.sh -q long
