# Path to local tools
export PATH=/groups/algorithm/tools/install/bin:$PATH

# Path can optionally be passed as an argument when auto detection fails
if [ "$#" -gt 0 ]
then
    config_dir="$1"
else
    config_dir=`dirname $BASH_ARGV[0]`
fi

# Figure out full path of config_dir
# Remove trailing / on config_dir because it will be added when used 
# by the script
if test ! -z "`uname | grep Darwin`"
then
    config_dir=`python -c 'import os,sys;print os.path.realpath(sys.argv[1])' $config_dir | sed 's/\/$//'`
else
    config_dir=`readlink -f $config_dir | sed 's/\/$//'`
fi

# Maybe you want to reference this directory in a script...
export L2_EXE_PATH=`dirname $config_dir`
export L2_FP_SRC_PATH=`dirname $config_dir`/lib/

# Remove this variable from environment
unset config_dir

# Set up environment based on hostname
curr_hostname=`hostname -s`

# Set up paths for OCO-2/ACOS machines
if [[ -e "/opt/local" ]]; then

    # List of packages to setup for use by L2 FP
    search_packages="valgrind gcc/4.8.1 gcc-4.5.1"

    for package in $search_packages; do
        pkg_base_dir=/opt/local/64/${package}

        if [ ! -e "$pkg_base_dir" ]; then
            pkg_base_dir=/opt/local/depot/${package}
        fi

        # Skip if package does not exist
        if [ ! -e "$pkg_base_dir" ]; then
            continue
        fi

        new_bin_path=${pkg_base_dir}/bin
        new_lib_path=${pkg_base_dir}/lib64

        if [ ! -e "$new_lib_path" ]; then
            new_lib_path=${pkg_base_dir}/lib
        fi

        for chk_path in $new_bin_path $new_lib_path; do
            if [ ! -e "$chk_path" ]; then
            echo "$chk_path does not exist!"
            exit 1
            fi
        done

        PATH=${new_bin_path}:${PATH}
        LD_LIBRARY_PATH=${new_lib_path}:${LD_LIBRARY_PATH}

        # man files not essential if not present
        new_man_path=${pkg_base_dir}/man
        if [ -e "$new_man_path" ]; then
            MANPATH=${new_man_path}:${MANPATH}
        fi

        unset pkg_base_dir new_bin_path new_man_path new_lib_path 
    done

    # Intel compiler has its own environmental variable script
    # Search first at scf-srv2 location then for scf-srv3 location
    ifort_base_dirs="/opt/local/depot/intel/Compiler/2015 /opt/local/depot/intel/Compiler/2013"
    for ifort_base in $ifort_base_dirs; do
        if [ -e "$ifort_base" ]; then
            # Set up first ifort compiler found
            source "$ifort_base/bin/ifortvars.sh" intel64
            break
        fi
    done
    
    unset search_packages

elif [[ "$curr_hostname" == "galaxy" || "$curr_hostname" == "nebula" ]]; then

    module load absoft-fortran-compiler.9.0 f_compiler/intel-ifort-compiler-10.1.018

fi # Other systems expect to have compilers set up system wide

# Remove these temporary variables
unset curr_hostname bits

export PATH LD_LIBRARY_PATH MANPATH
