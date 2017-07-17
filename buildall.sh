#!/bin/bash -e
export MROOT=`pwd`

cd core
make s p d r rs
cd ..

cd simp
make s p d r rs
cd ..