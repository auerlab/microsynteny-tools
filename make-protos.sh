#!/bin/sh -e

: ${PREFIX:=../local}
: ${LOCALBASE:=/usr/local}

# Undefined macro in stdlib.h
extras=-D_Noreturn=""

for file in *.c; do
    cproto $extras -I$LOCALBASE/include $file
done
