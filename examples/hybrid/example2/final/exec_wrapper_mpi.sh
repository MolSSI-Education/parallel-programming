#!/bin/bash
export MV2_ENABLE_AFFINITY=0
exec=$1
if [ $PMI_RANK -eq 0 ]; then
/usr/bin/time -v $exec
else
$exec
fi
