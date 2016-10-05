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

    # Fullerene has custom python location. Put before other paths
    # so this python version takes predence
    if [ -e "/opt/local/depot/python-64/2.7.3/bin" ]
    then
        export PATH="/opt/local/depot/python-64/2.7.3/bin:$PATH"
        export LD_LIBRARY_PATH="/opt/local/depot/python-64/2.7.3/lib:$LD_LIBRARY_PATH"
    fi

    # scf-srv3 has yet another python location
    if [ -e "/opt/local/depot/python/2.7.5/bin" ]
    then
        export PATH="/opt/local/depot/python/2.7.5/bin:$PATH"
        export LD_LIBRARY_PATH="/opt/local/depot/python/2.7.5/lib:$LD_LIBRARY_PATH"
    fi

    # Path to h5diff, the system one doesn't work for us
    if [ -e "/opt/local/depot/hdf5/1.8.7/bin" ]
    then
        export PATH="/opt/local/depot/hdf5/1.8.7/bin:$PATH"
        export LD_LIBRARY_PATH="/opt/local/depot/python/2.7.5/lib:$LD_LIBRARY_PATH"
    fi

    # Set up paths so executable scripts can be found
    export PATH="$PATH:$support_dir/utils"

    # Allow python modules to be found by utilities
    export PYTHONPATH="$PYTHONPATH:$support_dir"
fi

# Remove this variable from environment
unset support_dir
