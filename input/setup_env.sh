# Path can optionally be passed as an argument when auto detection fails
if [ "$#" -gt 0 ]
then
    input_dir="$1"
else
    input_dir=`dirname $BASH_ARGV[0]`
fi

# Figure out full path
# Remove trailing / because it will be added when used 
# by the script
if test ! -z "`uname | grep Darwin`"
then
    input_dir=`python -c 'import os,sys;print os.path.realpath(sys.argv[1])' $input_dir | sed 's/\/$//'`
else
    input_dir=`readlink -f $input_dir | sed 's/\/$//'`
fi

# Make sure that the directory is a valid one
if [ ! -e "$input_dir/setup_env.sh" ]
then
    echo "L2_Input: The input directory specified: '$input_dir' is not valid"
else
    export L2_INPUT_PATH="$input_dir"
fi

# Remove this variable from environment
unset input_dir
