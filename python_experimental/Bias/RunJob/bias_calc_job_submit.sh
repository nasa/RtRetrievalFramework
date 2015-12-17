#!/bin/bash
numsounding=`cat sounding_list.txt | wc -w`
endindex=`expr $numsounding - 1`
qsub -j oe -o log/bias_calc.log -t 0-${endindex} ./bias_calc_job.sh -q long
