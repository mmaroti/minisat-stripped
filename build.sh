#!/bin/bash -e
export MROOT=`pwd`

for DIR in core simp
do
  cd minisat/$DIR

  for BIN in s p d r
  do
    make $BIN
  done

  for BIN in debug profile release standard
  do
    make lib_$BIN.a
  done

  if [[ "$TRAVIS_OS_NAME" != "osx" ]]
  then
    make rs
  fi

  cd ../..
done
