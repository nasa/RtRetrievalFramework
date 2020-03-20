# Remove this variable from environment
if [ "$#" -gt 0 ]
then
    support_dir="$1"
else
    support_dir=`dirname $BASH_ARGV[0]`
fi

# Figure out full path
# Remove trailing / because it will be added when used 
# by the script
if test ! -z "`uname | grep Darwin`"
then
    support_dir=`python -c 'import os,sys;print os.path.realpath(sys.argv[1])' $support_dir | sed 's/\/$//'`
else
    support_dir=`readlink -f $support_dir | sed 's/\/$//'`
fi

# Make sure that the directory is a valid one
if [ ! -e "$support_dir/setup_env.sh" ]
then
    echo "L2_Support: The support directory specified: '$support_dir' is not valid"
else
    export L2_SUPPORT_PATH="$support_dir"
    export L2_SUPPORT_UTIL_PATH="$support_dir/utils"

    # Assume if a virtual environment is set that we want to use that
    # version of python. Otherwise, select the default python 3 version
    # We have a development version at /groups/algorithm/venv3, but
    # instead use the production version as the default
    DEFAULT_VENV="/groups/sdos/oco3/python_env/l2fp/20200303_l2fp"
    if ([ -z ${VIRTUAL_ENV+x} ] && [ -e "$DEFAULT_VENV/bin/activate" ])
    then
	VIRTUAL_ENV_DISABLE_PROMPT=t source $DEFAULT_VENV/bin/activate
    fi

    # Path to h5diff, the system one doesn't work for us
    if [ -e "/opt/local/depot/hdf5/1.8.14/bin" ]
    then
        export PATH="/opt/local/depot/hdf5/1.8.14/bin:$PATH"
    fi

    # Set up paths so executable scripts can be found
    export PATH="$PATH:$support_dir/utils"

    # Allow python modules to be found by utilities
    export PYTHONPATH="$PYTHONPATH:$support_dir"
fi

# Remove this variable from environment
unset support_dir
