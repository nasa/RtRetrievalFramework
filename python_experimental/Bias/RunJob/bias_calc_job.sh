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
sid="${id_list[${job_index}]}"
xco2_bias_calc --lua-config=config_quicklook_aerosol_test_0.lua --distribution-covariance=${sid}/cov_draw.pkl /groups/algorithm/l2_fp/oco2_simulator_test/automated/tags/B3.5.00_aerosol_testing_0/Orbit022/oco_Orbit022.config ${sid} ${sid}/bias_calc.mat




