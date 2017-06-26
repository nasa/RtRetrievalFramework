# Path can optionally be passed as an argument when auto detection fails
if [ "$#" -gt 0 ]
then
    package_dir="$1"
else
    package_dir=`dirname $BASH_ARGV[0]`
fi
if test ! -z "`uname | grep Darwin`"
then
    package_dir=`python -c 'import os,sys;print os.path.realpath(sys.argv[1])' $package_dir | sed 's/\/$//'`
else
    package_dir=`readlink -f $package_dir | sed 's/\/$//'`
fi

source $package_dir/config/setup_compiler_env.sh $package_dir/bin
source $package_dir/input/setup_env.sh $package_dir/input/
source $package_dir/support/setup_env.sh $package_dir/support
export L2_TCCON_SMALL_SET_PATH=$package_dir/test/tccon_small_set
export L2_TEST_PATH=$package_dir/test

if [ -n "$LUA_PATH" ]; then
    export LUA_PATH="$package_dir/input/common/config/?.lua;$package_dir/input/gosat/config/?.lua;$package_dir/input/oco/config/?.lua;$package_dir/input/fts/config/?.lua;${LUA_PATH}"
else
    export LUA_PATH="$package_dir/input/common/config/?.lua;$package_dir/input/gosat/config/?.lua;$package_dir/input/oco/config/?.lua;$package_dir/input/fts/config/?.lua;"
fi
# Remove this variable from environment
unset package_dir
