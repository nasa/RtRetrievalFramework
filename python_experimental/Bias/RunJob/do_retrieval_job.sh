#!/bin/bash
# Short script to run all the create_draw programs.
export PYTHONPATH=/home/smyth/Level2PythonBuild/install/lib/python2.7/site-packages:/home/smyth/Level2PythonBuild/install/lib/python2.7/site-packages:/home/smyth/Install/lib/python2.7/site-packages
export PATH=/home/smyth/Level2PythonBuild/install/bin:/groups/algorithm/tools/install/bin:/home/smyth/Install/bin:/opt/local/depot/python-64/2.7.3/bin:/opt/local/depot/gcc-4.5.1/bin:/opt/local/depot/emacs/23.2/bin:/opt/local/depot/intel/11.1/064/bin/intel64:/usr/NX/bin:/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/opt/local/64/bin:/opt/local/32/bin:/opt/dell/srvadmin/bin:/opt/local/64/bin:/home/smyth/bin
export LD_LIBRARY_PATH=/opt/local/depot/gcc-4.5.1/lib64:/opt/local/depot/python-64/2.7.3/lib:/groups/algorithm/tools/install/lib:/groups/algorithm/tools/install/lib64:/opt/local/depot/intel/11.1/064/lib/intel64:/opt/local/depot/intel/11.1/064/mkl/lib/em64t:/opt/local/64/lib
export LUA_PATH=/home/smyth/Level2PythonBuild/install/etc/full_physics/config/?.lua;
cd /scratch/algorithm/smyth/Bias
id_list=( `cat sounding_list.txt` )

if [[ ! -z "${PBS_ARRAYID}" ]]; then
   job_index=${PBS_ARRAYID}
else
   job_index=$1
fi
numdraw=700
sid_ind=`expr $job_index / $numdraw`
drawind=`expr $job_index % $numdraw`
sid="${id_list[${sid_ind}]}"
python ~/Level2/python_experimental/Bias/do_retrieval.py ${sid} ${sid}/cov_initial.shlv ${sid}/draw_${drawind}.shlv ${sid}/retrieval_initial_${drawind}.shlv

