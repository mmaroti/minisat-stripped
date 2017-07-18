#!/bin/bash -e
export MROOT=`pwd`

for DIR in core simp
do
  cd $DIR

  for BIN in p d r rs
  do
    make $BIN
  done

  for BIN in debug profile release standard
  do
    make lib_$BIN.a
  done

  if [[ "$TRAVIS_OS_NAME" != "osx" ]]
  then
    make s
  fi

  cd ..
done
