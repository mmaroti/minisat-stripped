#!/bin/bash -e
export MROOT=`pwd`

for DIR in core simp
do
  cd minisat/$DIR

  for BIN in s p d r rs libs libp libd libr
  do
    make $BIN
  done

  cd ../..
done
