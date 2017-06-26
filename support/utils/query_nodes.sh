#!/bin/sh



for node in `pbsnodes -a | grep -E '^\w+' | grep -v 'adm$'`
do
    echo $node
    printf "%"`echo $node | wc -c`"s\n" | sed 's| |-|g'

    ssh -o ForwardX11=no $node $*
    echo
done
