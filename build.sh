#!/bin/bash -e
export MROOT=`pwd`

for STD in c++98 c++11 c++14
do
  export CFLAGS="-std=$STD -Wall -Wno-parentheses -Wextra"
  
  if [ "$STD" != c++98 ]
  then
    export CFLAGS="$CFLAGS -Wmissing-noreturn"
  fi
  
  if [ "$CXX" == clang++ ]
  then
    export CFLAGS="$CFLAGS -Wweak-vtables"
  fi

  for DIR in core simp
  do
    echo =======================================
    echo ===== Building \'$DIR\' under $STD =====
    echo =======================================

    cd minisat/$DIR
    make clean
    rm -f *.a

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
done
