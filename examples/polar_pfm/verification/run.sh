#!/usr/bin/env bash

for var in upsilon theta
do
  for T in `seq 392 10 532`
  #for T in `seq 402 20 522`
  do
    ../../../magpie-opt -i velocity_${var}.i Temp_=${T} KE=1.933
  done
done
