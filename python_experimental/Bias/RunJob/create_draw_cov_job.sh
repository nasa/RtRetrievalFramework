#!/bin/bash
cd /scratch/algorithm/smyth/Bias
for sid in `cat sounding_list.txt`
do
    python ~/Level2/python_experimental/Bias/create_draw_cov.py ${sid} ${sid}/cov_draw.shlv ${sid}/cov_draw.pkl
done
