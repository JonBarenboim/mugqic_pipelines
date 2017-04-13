#!/bin/bash

# Setup environment to run MUGQIC pipelines
module load mugqic-pipelines/2.2.0
export PYTHONPATH=$MUGQIC_PIPELINES_HOME:$PYTHONPATH

$MUGQIC_PIPELINES_HOME/pipelines/episeq/episeq.py \
 -s 1-19 \
 -o output \
 -j pbs \
 -l debug \
 -d ./episeq.design \
 -r ./episeq.readset \
 -c episeq.ini 1> qsub.sh 2> debug.log

# Other flags
#--report || --clean

chmod +x qsub.sh

