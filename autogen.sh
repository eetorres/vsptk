#!/bin/bash

aclocal
autoheader
automake --add-missing
autoconf
#./configure CXX="icpc" CC="icc"
./configure
