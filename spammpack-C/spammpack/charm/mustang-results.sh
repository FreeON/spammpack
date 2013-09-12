#!/bin/bash

for i in slurm*out; do
  RESULT=$( grep -E "iteration [0-9]" $i | tail -n1 | awk '{print $4, $10}' )
  if [[ -z $RESULT ]]; then
    echo "no results in $i" 1>&2
  else
    echo $RESULT
  fi
done
