#!/bin/bash

# Read the functional (default to "b3lyp")
if [ ! -z "$1" ]; then
  functional="$1"
else
  read -p 'Functional? (Defaults to b3lyp)' functional
fi

if [ -z $functional ]; then
  functional="b3lyp"
fi

for molecule_name in $(cat < molecule_name.txt)
do  
    #$1: functional
    python generate_gjf_from_name.py $molecule_name $functional &

    for job in $(jobs -p)
    do
        wait $job
    done
done

