#!/bin/bash

# Wrapper to pbsdsh python scripts. Exports
# required env variables.

# set env variables

PATH=$1:PATH
export LD_LIBRARY_PATH=$2

# set python arguments
exe=$3
output=$4
classname=$5
funcname=$6
pythonpath=$7

$exe $output $classname $funcname $pythonpath
