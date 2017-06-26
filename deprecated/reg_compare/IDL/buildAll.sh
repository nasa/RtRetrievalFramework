#!/bin/sh

echo
echo "OCO L2 Matrix Compare"
echo "--------------"
cd OCOL2MatrixCompare
echo '@make_idl' | idl -32
cd ..

