#!/bin/bash

function Usage
{
    echo "Usage: $0 -p project -r rev1 -rev2 -t tag1 -t tag2 [-f file]"
    echo "either -r or -t can be specified"
    echo "omitting -f will create a list of differing files"
    echo "e.g. $0 -p L2_EXE -t L2_20060307 -t L2_20060411 -f src/modules/ainput.F90"
    echo "     $0 -p L2_EXE -r 36 -r 38"
    exit 1
}

function svn_rev_latest
{
    project=$1
    rev=$(svn info $SVNROOT/$project | grep Rev: | cut -d':' -f 2 | sed -e 's/ //')
}

function svn_rev_tag
{
    project=$1
    tag=$2
    SVN_REV_TAG_VALUE=$(svn info $SVNROOT/$project/tags/$tag | grep Rev: | cut -d':' -f 2 | awk '{print $1}')
}

# No arguments specified, so show usage
if [ -z $1 ]; then
    Usage
fi

if [ -z "$SVNROOT" ]; then
    echo "\$SVNROOT must be defined!"
    echo "e.g. export SVNROOT=https://svn/oco/alg/"
    exit 1
fi

# parse arguments
while getopts "d:f:p:r:t:" opt; do
    case $opt in
	d ) dir=$OPTARG
	    subdir=$OPTARG ;;
	f ) file=$OPTARG ;;
	p ) project=$OPTARG ;;
	r ) if [ -z $rev1 ]; then
	    rev1=$OPTARG
	    else
	    rev2=$OPTARG
	    fi ;;
	t ) if [ -z $tag1 ]; then
	    tag1=$OPTARG
	    else
	    tag2=$OPTARG
	    fi ;;
	* ) Usage ;;
    esac
done

if [ -z $subdir ]; then
    dir="trunk"
else
    dir="trunk/$subdir"
fi

# Get the two revision numbers to compare
# If rev1 is not defined, get it from tag1
if [ -z $rev1 ]; then
    if [ -z $tag1 ]; then
	echo "error: neither -r nor -t was used"
	exit 1
    fi
    svn_rev_tag $project $tag1
    rev1=$SVN_REV_TAG_VALUE
fi

# if rev2 is not defined, get it from tag1.  If tag1 was used for rev1
# then use tag2.  If neither tag1 nor tag2 was not specified, use the
# latest revision.
if [ -z $rev2 ]; then
    if [ -z $tag1 ]; then
	# Use the latest revision
	svn_rev_latest $project
	rev2=$SVN_REV_TAG_VALUE
    else # tag1 was specified
	svn_rev_tag $project $tag1
	rev2=$SVN_REV_TAG_VALUE
	if [ $rev2 -eq $rev1 ]; then
	    if [ -z $tag2 ]; then
	        # Use the latest revision
		svn_rev_latest $project
		rev2=$SVN_REV_TAG_VALUE
	    else
		svn_rev_tag $project $tag2
		rev2=$SVN_REV_TAG_VALUE
	    fi
	fi
    fi
fi

# If no file was specified, create a list of changed files
if [ -z $file ]; then
    string="tkdiffrevs.sh -p $project -r $rev1 -r $rev2 -f $subdir/"
    svn diff $SVNROOT/$project/$dir@$rev1 $SVNROOT/$project/$dir@$rev2 | grep Index: | sort | sed -e "s,Index: ,$string," 
    exit 1
fi

echo svn diff --diff-cmd tkdiffwrap.sh $SVNROOT/$project/$dir/$file@$rev1 $SVNROOT/$project/$dir/$file@$rev2

svn diff --diff-cmd tkdiffwrap.sh $SVNROOT/$project/$dir/$file@$rev1 $SVNROOT/$project/$dir/$file@$rev2
